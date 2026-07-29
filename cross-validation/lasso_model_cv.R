# =============================================================================
# Norovirus LASSO (glmnet) - Spatial Cross-Validation
# =============================================================================

setwd("~/GitHub/Statistical-Learning-Methods-for-Wastewater-Pathogen-Prediction")


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(glmnet)
library(sf)
library(sp)
library(raster)
library(blockCV)

set.seed(123)

# -----------------------------------------------------------------------------
# Load and Prepare Data
# -----------------------------------------------------------------------------

nov_df <- readr::read_csv("data/nov_df.csv") %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
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
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700) %>%
  st_transform(crs = bng_km) %>%
  arrange(site_code) %>%
  mutate(
    site_code = as.factor(site_code),
    s_index   = as.numeric(site_code),
    site_code = as.character(site_code)
  )

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
# Model Variables
# -----------------------------------------------------------------------------

vars_used <- c(
  "avg_Log10_NoV", "lockdown_step3", "lockdown_step4", "lockdown_planB",
  "scale_school_den", "scale_carehome_den", "scale_mob_7day",
  "scale_bame", "scale_imd", "scale_prop_agri", "scale_prop_urb",
  "temp_2", "temp_3", "temp_4", "temp_5",
  "rain_2", "rain_3", "rain_4", "rain_5",
  "week", "site_code"
)

model_rhs <- ~ -1 +
  lockdown_step3 + lockdown_step4 + lockdown_planB +
  scale_school_den + scale_carehome_den + scale_mob_7day +
  scale_bame + scale_imd + scale_prop_agri + scale_prop_urb +
  temp_2 + temp_3 + temp_4 + temp_5 +
  rain_2 + rain_3 + rain_4 + rain_5 +
  week * site_code

# -----------------------------------------------------------------------------
# Spatial Cross-Validation Loop
# -----------------------------------------------------------------------------

fit         <- list()
metrics     <- data.frame()
predictions <- data.frame()

# Start timer
cv_start_time <- proc.time()

for (k in 1:10) {
  
  # Split into training and validation sets
  train <- subset(nov, folds != k) %>%
    st_drop_geometry() %>%
    dplyr::filter(complete.cases(dplyr::select(., all_of(vars_used))))
  
  val <- subset(nov, folds == k) %>%
    st_drop_geometry() %>%
    dplyr::filter(complete.cases(dplyr::select(., all_of(vars_used))))
  
  # Build model matrices
  y_train <- train$avg_Log10_NoV
  y_test  <- val$avg_Log10_NoV
  X_train <- model.matrix(model_rhs, data = train)
  X_test  <- model.matrix(model_rhs, data = val)
  
  # Align test matrix columns to training matrix
  missing_cols <- setdiff(colnames(X_train), colnames(X_test))
  extra_cols   <- setdiff(colnames(X_test),  colnames(X_train))
  
  if (length(extra_cols)   > 0) X_test <- X_test[, !colnames(X_test) %in% extra_cols, drop = FALSE]
  if (length(missing_cols) > 0) {
    zero_mat <- matrix(0, nrow = nrow(X_test), ncol = length(missing_cols),
                       dimnames = list(NULL, missing_cols))
    X_test <- cbind(X_test, zero_mat)
  }
  X_test <- X_test[, colnames(X_train), drop = FALSE]
  
  # Fit LASSO with cross-validated lambda
  fit[[k]]    <- cv.glmnet(X_train, y_train, alpha = 1, family = "gaussian")
  best_lambda <- fit[[k]]$lambda.min
  
  val$predicted <- as.numeric(predict(fit[[k]], s = best_lambda, newx = X_test))
  
  # ---------------------------------------------------------------------------
  # Bootstrap Uncertainty Estimation (n = 500)
  # ---------------------------------------------------------------------------
  
  n_boot     <- 500
  boot_preds <- matrix(NA, nrow = n_boot, ncol = nrow(val))
  
  for (i in 1:n_boot) {
    set.seed(123 + i)
    boot_idx       <- sample(seq_len(nrow(train)), replace = TRUE)
    boot_model     <- glmnet(X_train[boot_idx, ], y_train[boot_idx],
                             alpha = 1, lambda = best_lambda, family = "gaussian")
    boot_preds[i, ] <- as.numeric(predict(boot_model, s = best_lambda, newx = X_test))
  }
  
  val$q0.025 <- apply(boot_preds, 2, quantile, probs = 0.025)
  val$q0.975 <- apply(boot_preds, 2, quantile, probs = 0.975)
  
  val_out <- val %>%
    dplyr::select(site_code, s_index, week_beginning, week,
                  avg_Log10_NoV, predicted, q0.025, q0.975)
  
  predictions <- rbind(predictions, val_out)
  
  # ---------------------------------------------------------------------------
  # Compute Fold-Level Performance Metrics
  # ---------------------------------------------------------------------------
  
  resid <- val$avg_Log10_NoV - val$predicted
  mse   <- mean(resid^2,           na.rm = TRUE)
  rmse  <- sqrt(mse)
  mae   <- mean(abs(resid),         na.rm = TRUE)
  mape  <- mean(abs(resid / val$avg_Log10_NoV), na.rm = TRUE) * 100
  bias  <- mean(-resid,             na.rm = TRUE)
  pbias <- (bias / mean(val$avg_Log10_NoV, na.rm = TRUE)) * 100
  corr  <- as.numeric(cor(val$avg_Log10_NoV, val$predicted,
                          use = "complete.obs", method = "pearson"))
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
  
  metrics <- rbind(metrics, metrics_fold)
  message(paste(k, "folds done"))
}

cv_elapsed <- proc.time() - cv_start_time

elapsed_sec <- cv_elapsed["elapsed"]
hours       <- floor(elapsed_sec / 3600)
minutes     <- floor((elapsed_sec %% 3600) / 60)
seconds     <- round(elapsed_sec %% 60, 1)

message(sprintf("CV loop completed in %dh %dm %.1fs", hours, minutes, seconds))

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

