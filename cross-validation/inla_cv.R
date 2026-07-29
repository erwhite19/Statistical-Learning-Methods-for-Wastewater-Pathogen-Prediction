# =============================================================================
# Norovirus INLA Spatio-Temporal Model - Spatial Cross-Validation
# =============================================================================

setwd("~/GitHub/Statistical-Learning-Methods-for-Wastewater-Pathogen-Prediction")

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(INLA)
library(inlabru)
library(sf)
library(sp)
library(blockCV)
library(raster)

set.seed(123)

# -----------------------------------------------------------------------------
# Load and Prepare Data
# -----------------------------------------------------------------------------

nov_df <- readr::read_csv("data/nov_df.csv") %>%
  dplyr::mutate(row_id = dplyr::row_number())

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
COV <- function(z, lower = NULL, upper = NULL) {
  x <- as.matrix(z) >= as.matrix(lower) & as.matrix(z) <= as.matrix(upper)
  u <- x[!is.na(x)]
  round(sum(u) / length(u) * 100, 4)
}

# Z-score standardisation applied to specified covariates
scale_covariates <- function(df, covariates_to_scale) {
  scale_fn <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  df %>% mutate(across(all_of(covariates_to_scale), scale_fn, .names = "scale_{.col}"))
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
    site_code = as.character(site_code),
    Week_date = as.Date(week_beginning)
  ) %>%
  arrange(Week_date) %>%
  mutate(f_index = as.numeric(Week_date))

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
# INLA Priors and Model Formula
# -----------------------------------------------------------------------------

pc_prec <- list(prec = list(prior = "pc.prec", param = c(10, 0.05)))

formula <- as.formula(
  'avg_Log10_NoV ~ -1 + Intercept +
   lockdown_step3 + lockdown_step4 + lockdown_planB +
   scale_school_den + scale_carehome_den + scale_mob_7day +
   scale_bame + scale_imd + scale_prop_agri + scale_prop_urb +
   temp_2 + temp_3 + temp_4 + temp_5 +
   rain_2 + rain_3 + rain_4 + rain_5 +
   f(site_code, model = "iid", hyper = pc_prec) +
   f(week,      model = "iid", hyper = pc_prec) +
   f(spatial.field, model = spde,
     group         = spatial.field.group,
     control.group = list(model = "iid", hyper = pc_prec))'
)

# -----------------------------------------------------------------------------
# Spatial Cross-Validation Loop
# -----------------------------------------------------------------------------

metrics     <- data.frame()
predictions <- data.frame()

# Start timer
cv_start_time <- proc.time()

for (k in 1:10) {
  
  cat("\n========== Fold", k, "of 10 ==========\n")
  
  # Split into training and validation sets (sf for mesh, data frame for stacks)
  train_sf <- nov %>% filter(folds != k)
  val_sf   <- nov %>% filter(folds == k)
  train    <- st_drop_geometry(train_sf)
  val      <- st_drop_geometry(val_sf)
  
  # ---------------------------------------------------------------------------
  # Build Mesh on Training Data Only
  # ---------------------------------------------------------------------------
  
  coords_train <- as.matrix(st_coordinates(train_sf))
  max.edge     <- 10
  bound.outer  <- diff(range(coords_train[, 1])) / 5
  domain       <- inla.nonconvex.hull(coords_train)
  
  mesh <- fm_mesh_2d_inla(
    boundary = domain,
    loc      = coords_train,
    max.edge = c(1, 2) * max.edge,
    offset   = c(max.edge, bound.outer),
    cutoff   = max.edge / 5,
    crs      = st_crs(train_sf)
  )
  
  # ---------------------------------------------------------------------------
  # SPDE Object
  # ---------------------------------------------------------------------------
  
  spde <- inla.spde2.pcmatern(
    mesh        = mesh,
    alpha       = 2,
    prior.range = c(10, 0.05),
    prior.sigma = c(5,  0.05)
  )
  
  # ---------------------------------------------------------------------------
  # Estimation Stack (Training Data)
  # ---------------------------------------------------------------------------
  
  n_week_train <- length(unique(train$week))
  
  s.index_est <- inla.spde.make.index(
    name    = "spatial.field",
    group   = train$week,
    n.spde  = spde$n.spde,
    n.group = n_week_train
  )
  
  A_est <- inla.spde.make.A(
    mesh    = mesh,
    loc     = coords_train,
    group   = train$week,
    n.group = n_week_train
  )
  
  stack_est <- inla.stack(
    tag  = "est",
    data = list(avg_Log10_NoV = train$avg_Log10_NoV),
    A    = list(A_est, 1),
    effects = list(
      c(s.index_est, list(Intercept = 1)),
      list(
        week               = train$week,
        site_code          = train$site_code,
        lockdown_step3     = train$lockdown_step3,
        lockdown_step4     = train$lockdown_step4,
        lockdown_planB     = train$lockdown_planB,
        lockdown_lifting   = train$lockdown_lifting,
        scale_school_den   = train$scale_school_den,
        scale_carehome_den = train$scale_carehome_den,
        scale_mob_7day     = train$scale_mob_7day,
        scale_bame         = train$scale_bame,
        scale_imd          = train$scale_imd,
        scale_prop_urb     = train$scale_prop_urb,
        scale_prop_agri    = train$scale_prop_agri,
        temp_2 = train$temp_2, temp_3 = train$temp_3,
        temp_4 = train$temp_4, temp_5 = train$temp_5,
        rain_2 = train$rain_2, rain_3 = train$rain_3,
        rain_4 = train$rain_4, rain_5 = train$rain_5
      )
    )
  )
  
  # ---------------------------------------------------------------------------
  # Validation Stack (Response Set to NA for Prediction)
  # ---------------------------------------------------------------------------
  
  coords_val  <- as.matrix(st_coordinates(val_sf))
  n_week_val  <- length(unique(val$week))
  
  s.index_val <- inla.spde.make.index(
    name    = "spatial.field",
    n.spde  = spde$n.spde,
    n.group = n_week_val
  )
  
  A_val <- inla.spde.make.A(
    mesh    = mesh,
    loc     = coords_val,
    group   = val$week,
    n.group = n_week_val
  )
  
  stack_val <- inla.stack(
    tag  = "val",
    data = list(avg_Log10_NoV = NA),
    A    = list(A_val, 1),
    effects = list(
      c(s.index_val, list(Intercept = 1)),
      list(
        week               = val$week,
        site_code          = val$site_code,
        lockdown_step3     = val$lockdown_step3,
        lockdown_step4     = val$lockdown_step4,
        lockdown_planB     = val$lockdown_planB,
        lockdown_lifting   = val$lockdown_lifting,
        scale_school_den   = val$scale_school_den,
        scale_carehome_den = val$scale_carehome_den,
        scale_mob_7day     = val$scale_mob_7day,
        scale_bame         = val$scale_bame,
        scale_imd          = val$scale_imd,
        scale_prop_urb     = val$scale_prop_urb,
        scale_prop_agri    = val$scale_prop_agri,
        temp_2 = val$temp_2, temp_3 = val$temp_3,
        temp_4 = val$temp_4, temp_5 = val$temp_5,
        rain_2 = val$rain_2, rain_3 = val$rain_3,
        rain_4 = val$rain_4, rain_5 = val$rain_5
      )
    )
  )
  
  # ---------------------------------------------------------------------------
  # Join Stacks and Fit INLA Model
  # ---------------------------------------------------------------------------
  
  join_stack <- inla.stack(stack_est, stack_val)
  
  fit_inla <- inla(
    formula,
    data              = inla.stack.data(join_stack),
    family            = "gaussian",
    control.predictor = list(
      A       = inla.stack.A(join_stack),
      compute = TRUE,
      link    = 1
    ),
    control.compute = list(
      dic                        = TRUE,
      waic                       = TRUE,
      cpo                        = TRUE,
      config                     = TRUE,
      return.marginals.predictor = TRUE
    ),
    control.fixed = list(
      mean           = 0,
      prec           = 0.0001,
      mean.intercept = 0,
      prec.intercept = 0.0001
    ),
    control.family = list(hyper = list(prec = list(param = c(1, 0.2)))),
    verbose = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Extract Validation Predictions
  # ---------------------------------------------------------------------------
  
  val_idx       <- inla.stack.index(join_stack, tag = "val")$data
  val$predicted <- fit_inla$summary.fitted.values[val_idx, "mean"]
  val$sd        <- fit_inla$summary.fitted.values[val_idx, "sd"]
  val$q0.025    <- fit_inla$summary.fitted.values[val_idx, "0.025quant"]
  val$q0.975    <- fit_inla$summary.fitted.values[val_idx, "0.975quant"]
  
  val_out <- val %>%
    dplyr::select(site_code, week, avg_Log10_NoV, predicted, sd, q0.025, q0.975)
  
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
  cat("Fold", k, "complete | RMSE:", round(rmse, 4),
      "| Coverage:", round(cov, 2), "%\n")
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

cat("Overall Coverage Probability (%):", coverage_probability, "\n")

# -----------------------------------------------------------------------------
# Summary Metrics Across Folds
# -----------------------------------------------------------------------------

summary_metrics <- metrics %>%
  summarise(
    Mean_MSE      = mean(MSE),
    SE_MSE        = sd(MSE)   / sqrt(n()),
    CI_Lower_MSE  = Mean_MSE  - qt(0.975, df = n() - 1) * SE_MSE,
    CI_Upper_MSE  = Mean_MSE  + qt(0.975, df = n() - 1) * SE_MSE,
    Mean_RMSE     = mean(RMSE),
    SE_RMSE       = sd(RMSE)  / sqrt(n()),
    CI_Lower_RMSE = Mean_RMSE - qt(0.975, df = n() - 1) * SE_RMSE,
    CI_Upper_RMSE = Mean_RMSE + qt(0.975, df = n() - 1) * SE_RMSE,
    Mean_MAE      = mean(MAE),
    Mean_MAPE     = mean(MAPE),
    Mean_BIAS     = mean(BIAS),
    Mean_pBIAS    = mean(pBIAS),
    Mean_CORR     = mean(CORR),
    COV           = coverage_probability
  )

print(summary_metrics)

