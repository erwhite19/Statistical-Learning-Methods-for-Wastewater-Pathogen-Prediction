# =============================================================================
# Norovirus GAM (mgcv) - Spatial Cross-Validation
# =============================================================================

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(mgcv)
library(sf)
library(sp)
library(raster)
library(blockCV)

set.seed(123)

# -----------------------------------------------------------------------------
# Load and Prepare Data
# -----------------------------------------------------------------------------

nov_df <- read.csv("data/nov_df.csv") %>%
  dplyr::filter(!is.na(avg_Log10_NoV))

# Ensure no missing coordinates
coords_master <- nov_df %>%
  filter(!is.na(Easting) & !is.na(Northing)) %>%
  distinct(site_code, Easting, Northing)

nov_df <- nov_df %>%
  dplyr::select(-Easting, -Northing) %>%
  left_join(coords_master, by = "site_code")

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Coverage probability: proportion of observed values within predicted interval
COV <- function(z, lower = NULL, upper = NULL, coverage = NULL) {
  if (!is.null(lower) && !is.null(upper)) {
    x <- as.matrix(z) >= as.matrix(lower) & as.matrix(z) <= as.matrix(upper)
    u <- x[!is.na(x)]
    round(sum(u) / length(u) * 100, 4)
  } else if (!is.null(coverage)) {
    round(mean(coverage, na.rm = TRUE), 4)
  }
}

# Z-score standardisation applied to specified covariates
scale_covariates <- function(df, covariates_to_scale) {
  scale <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  df %>% mutate(across(all_of(covariates_to_scale), scale, .names = "scale_{.col}"))
}

# -----------------------------------------------------------------------------
# Scale Covariates and Project to km CRS
# -----------------------------------------------------------------------------

covariates_to_scale <- c(
  "school_den", "carehome_den", "imd", "bame",
  "mob_7day", "prop_agri", "prop_urb"
)

bng_km <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

nov <- scale_covariates(nov_df, covariates_to_scale) %>%
  as.data.frame() %>%
  arrange(site_code) %>%
  mutate(
    site_code = as.factor(site_code),
    s_index   = as.numeric(site_code),
    site_code = as.character(site_code),
    Week_date = as.Date(week_beginning)
  ) %>%
  arrange(Week_date) %>%
  mutate(f_index = as.numeric(Week_date)) %>%
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700) %>%
  st_transform(crs = bng_km)

# -----------------------------------------------------------------------------
# Build Spatial Fold Blocks (blockCV)
# -----------------------------------------------------------------------------

# Unique monitoring sites for fold assignment
sites <- nov %>%
  dplyr::select(site_code) %>%
  distinct()

england <- st_read("data/england_crop.shp") %>%
  st_transform(crs = st_crs(nov))

england_raster <- raster(england)

folds <- cv_spatial(
  x    = sites,
  r    = england_raster,
  k    = 10,
  seed = 12,
  plot = FALSE
)

fold_blocks <- folds$blocks %>%
  st_transform(crs = st_crs(nov))

# Assign fold membership to observations via spatial intersection
nov <- st_intersection(nov, fold_blocks)

# -----------------------------------------------------------------------------
# Spatial Cross-Validation Loop
# -----------------------------------------------------------------------------

fit         <- list()
metrics     <- data.frame()
predictions <- data.frame()

for (k in 1:10) {
  
  # Split into training and validation sets
  train <- subset(nov, folds != k)
  val   <- subset(nov, folds == k)
  
  # Extract coordinates back from geometry
  train_coords   <- st_coordinates(train)
  train$Easting  <- train_coords[, 1]
  train$Northing <- train_coords[, 2]
  
  val_coords   <- st_coordinates(val)
  val$Easting  <- val_coords[, 1]
  val$Northing <- val_coords[, 2]
  
  # Fit GAM with spatial and temporal smooths
  fit[[k]] <- gam(
    avg_Log10_NoV ~ -1 +
      s(Easting, Northing) + s(week) +
      lockdown_step3 + lockdown_step4 + lockdown_planB +
      scale_school_den + scale_carehome_den + scale_mob_7day +
      scale_bame + scale_imd + scale_prop_agri + scale_prop_urb +
      temp_2 + temp_3 + temp_4 + temp_5 +
      rain_2 + rain_3 + rain_4 + rain_5,
    data   = train,
    family = gaussian(),
    method = "REML"
  )
  
  print(summary(fit[[k]]))
  print(gam.check(fit[[k]]))
  
  # Predict on validation set with 95% confidence interval
  pred       <- predict(fit[[k]], newdata = val, type = "response", se.fit = TRUE)
  val$predicted <- pred$fit
  val$se.fit    <- pred$se.fit
  val$q0.025    <- pred$fit - 1.96 * pred$se.fit
  val$q0.975    <- pred$fit + 1.96 * pred$se.fit
  
  val <- val %>%
    dplyr::select(site_code, s_index, week, Week_date, avg_Log10_NoV, predicted, q0.025, q0.975) %>%
    st_drop_geometry()
  
  predictions <- rbind(predictions, val)
  
  # Compute fold-level performance metrics
  resid <- val$avg_Log10_NoV - val$predicted
  mse   <- mean(resid^2,           na.rm = TRUE)
  rmse  <- sqrt(mse)
  mae   <- mean(abs(resid),         na.rm = TRUE)
  mape  <- mean(abs(resid / val$avg_Log10_NoV), na.rm = TRUE) * 100
  bias  <- mean(-resid,             na.rm = TRUE)
  pbias <- (bias / mean(val$avg_Log10_NoV, na.rm = TRUE)) * 100
  corr  <- cor(val$avg_Log10_NoV, val$predicted, use = "complete.obs", method = "pearson")
  cov   <- COV(val$avg_Log10_NoV, lower = val$q0.025, upper = val$q0.975)
  
  metrics_fold <- data.frame(
    Fold  = k,
    MSE   = mse,
    RMSE  = rmse,
    MAE   = mae,
    MAPE  = mape,
    BIAS  = bias,
    pBIAS = pbias,
    CORR  = corr,
    COV   = cov
  )
  
  print(metrics_fold)
  metrics <- rbind(metrics, metrics_fold)
  message(paste(k, "folds done"))
}

# -----------------------------------------------------------------------------
# Overall Coverage Probability (Vectorised)
# -----------------------------------------------------------------------------

coverage_probability <- mean(
  predictions$avg_Log10_NoV >= predictions$q0.025 &
    predictions$avg_Log10_NoV <= predictions$q0.975,
  na.rm = TRUE
) * 100

print(paste("Coverage Probability (%):", coverage_probability))

# -----------------------------------------------------------------------------
# Summary Metrics Across Folds
# -----------------------------------------------------------------------------

summary_metrics <- metrics %>%
  summarise(
    Mean_MSE      = mean(MSE,  na.rm = TRUE),
    SE_MSE        = sd(MSE,    na.rm = TRUE) / sqrt(n()),
    CI_Lower_MSE  = Mean_MSE  - qt(0.975, df = n() - 1) * SE_MSE,
    CI_Upper_MSE  = Mean_MSE  + qt(0.975, df = n() - 1) * SE_MSE,
    Mean_RMSE     = mean(RMSE, na.rm = TRUE),
    SE_RMSE       = sd(RMSE,   na.rm = TRUE) / sqrt(n()),
    CI_Lower_RMSE = Mean_RMSE - qt(0.975, df = n() - 1) * SE_RMSE,
    CI_Upper_RMSE = Mean_RMSE + qt(0.975, df = n() - 1) * SE_RMSE,
    Mean_MAE      = mean(MAE,  na.rm = TRUE),
    Mean_MAPE     = mean(MAPE, na.rm = TRUE),
    Mean_BIAS     = mean(BIAS, na.rm = TRUE),
    Mean_pBIAS    = mean(pBIAS, na.rm = TRUE),
    Mean_CORR     = mean(CORR, na.rm = TRUE),
    COV           = coverage_probability
  )

print(summary_metrics)

# -----------------------------------------------------------------------------
# Save Outputs
# -----------------------------------------------------------------------------

write.csv(metrics,     "outputs/gam_cv_metrics_by_fold.csv",      row.names = FALSE)
write.csv(predictions, "outputs/gam_cv_predictions_by_fold.csv",  row.names = FALSE)