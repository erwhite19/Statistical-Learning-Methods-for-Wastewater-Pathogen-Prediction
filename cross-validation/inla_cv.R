library(tidyverse)
library(INLA)
library(inlabru)
library(sf)
library(sp)
library(blockCV)
library(raster)

set.seed(123)

# ── Load and prepare data (identical to your RF script) ──────────────────────

nov_df <- readr::read_csv("data/nov_df.csv") %>%
  dplyr::mutate(row_id = dplyr::row_number())
coords_master <- nov_df %>%
  filter(!is.na(Easting) & !is.na(Northing)) %>%
  distinct(site_code, Easting, Northing)

nov_df <- nov_df %>%
  dplyr::select(-Easting, -Northing) %>%
  left_join(coords_master, by = "site_code")

# Coverage function
COV <- function(z, lower = NULL, upper = NULL) {
  z     <- as.matrix(z)
  lower <- as.matrix(lower)
  upper <- as.matrix(upper)
  x <- z >= lower & z <= upper
  u <- x[!is.na(x)]
  round(sum(u) / length(u) * 100, 4)
}

# Scale covariates
covariates_to_scale <- c("school_den", "carehome_den", "imd", "bame",
                         "mob_7day", "prop_agri", "prop_urb")
scale_covariates <- function(df, covariates_to_scale) {
  scale_fn <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  df %>% mutate(across(all_of(covariates_to_scale), scale_fn,
                       .names = "scale_{.col}"))
}

nov_df <- scale_covariates(nov_df, covariates_to_scale)
nov_df <- as.data.frame(nov_df)

nov <- nov_df
nov <- st_as_sf(nov_df, coords = c("Easting", "Northing"), crs = 27700)
nov <- st_transform(nov, crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717
                    +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")

nov <- nov %>%
  arrange(site_code) %>%
  mutate(site_code = as.factor(site_code),
         s_index   = as.numeric(site_code),
         site_code = as.character(site_code))

nov$Week_date <- as.Date(nov$week_beginning)
nov <- nov %>%
  arrange(Week_date) %>%
  mutate(f_index = as.numeric(Week_date))



# ── Spatial folds (identical to RF) ──────────────────────────────────────────

sites <- nov %>% dplyr::select(site_code) %>% unique()
england <- st_read("data/england_crop.shp")
england <- st_transform(england, crs = st_crs(nov))
raster  <- raster(england)

folds <- cv_spatial(x = sites, r = raster, k = 10, seed = 12, plot = FALSE)
fold_blocks <- folds$blocks
fold_blocks  <- st_transform(fold_blocks, crs = st_crs(nov))
nov <- nov %>% st_intersection(fold_blocks)

# ── INLA priors and formula (identical to your full model) ───────────────────

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
     group = spatial.field.group,
     control.group = list(model = "iid", hyper = pc_prec))'
)

vars_used <- c(
  "avg_Log10_NoV", "lockdown_step3", "lockdown_step4", "lockdown_planB",
  "scale_school_den", "scale_carehome_den", "scale_mob_7day",
  "scale_bame", "scale_imd", "scale_prop_agri", "scale_prop_urb",
  "temp_2", "temp_3", "temp_4", "temp_5",
  "rain_2", "rain_3", "rain_4", "rain_5", "week", "site_code"
)

# ── CV loop ───────────────────────────────────────────────────────────────────

metrics     <- data.frame()
predictions <- data.frame()




for (k in 1:10) {
  
  cat("\n========== Fold", k, "of 10 ==========\n")
  
  # ── Split ──────────────────────────────────────────────────────────────────
  train <- subset(nov, folds != k) %>%
    st_drop_geometry() 
  
  val <- subset(nov, folds == k) %>%
    st_drop_geometry() 
  
  # Need coordinates back for mesh
  train_sf <- nov %>% filter(folds != k) 
  
  val_sf <- nov %>% filter(folds == k) 
  
  # ── Mesh (built on TRAINING data only) ────────────────────────────────────
  # Important: mesh built on training coords only — val sites are unobserved
  coords_train <- as.matrix(st_coordinates(train_sf))
  
  max.edge    <- 10
  bound.outer <- diff(range(coords_train[, 1])) / 5
  domain      <- inla.nonconvex.hull(coords_train)
  
  mesh <- fm_mesh_2d_inla(
    boundary = domain,
    loc      = coords_train,
    max.edge = c(1, 2) * max.edge,
    offset   = c(max.edge, bound.outer),
    cutoff   = max.edge / 5,
    crs      = st_crs(train_sf)
  )
  
  # ── SPDE ───────────────────────────────────────────────────────────────────
  spde <- inla.spde2.pcmatern(
    mesh        = mesh,
    alpha       = 2,
    prior.range = c(10, 0.05),
    prior.sigma = c(5,  0.05)
  )
  
  # ── Estimation stack (training data) ──────────────────────────────────────
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
  
  # ── Validation stack (val sites — response set to NA) ─────────────────────
  coords_val   <- as.matrix(st_coordinates(val_sf)) 
  n_week_val   <- length(unique(val$week))
  
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
    data = list(avg_Log10_NoV = NA),   # ← NA so INLA predicts here
    A    = list(A_val, 1),
    effects = list(
      c(s.index_val, list(Intercept = 1)),
      list(
        week               = val$week,
        site_code          = val$site_code,  # unseen sites → iid = 0
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
  
  # ── Join stacks and fit ────────────────────────────────────────────────────
  join_stack <- inla.stack(stack_est, stack_val)
  
  fit_inla <- inla(
    formula,
    data    = inla.stack.data(join_stack),
    family  = "gaussian",
    control.predictor = list(
      A       = inla.stack.A(join_stack),
      compute = TRUE,
      link    = 1
    ),
    control.compute = list(
      dic      = TRUE,
      waic     = TRUE,
      cpo      = TRUE,
      config   = TRUE,
      return.marginals.predictor = TRUE
    ),
    control.fixed  = list(mean = 0, prec = 0.0001,
                          mean.intercept = 0, prec.intercept = 0.0001),
    control.family = list(hyper = list(
      prec = list(param = c(1, 0.2))
    )),
    verbose = FALSE
  )
  
  # ── Extract validation predictions ────────────────────────────────────────
  # Get indices for the "val" tag in the joined stack
  val_idx <- inla.stack.index(join_stack, tag = "val")$data
  
  val$predicted <- fit_inla$summary.fitted.values[val_idx, "mean"]
  val$sd        <- fit_inla$summary.fitted.values[val_idx, "sd"]
  val$q0.025    <- fit_inla$summary.fitted.values[val_idx, "0.025quant"]
  val$q0.975    <- fit_inla$summary.fitted.values[val_idx, "0.975quant"]
  
  # ── Save predictions ───────────────────────────────────────────────────────
  val_out <- val %>%
    dplyr::select(site_code, week, avg_Log10_NoV,
                  predicted, sd, q0.025, q0.975)
  
  predictions <- rbind(predictions, val_out)
  
  # ── Metrics ────────────────────────────────────────────────────────────────
  mse   <- mean((val$avg_Log10_NoV - val$predicted)^2, na.rm = TRUE)
  rmse  <- sqrt(mse)
  mae   <- mean(abs(val$avg_Log10_NoV - val$predicted), na.rm = TRUE)
  mape  <- mean(abs((val$avg_Log10_NoV - val$predicted) /
                      val$avg_Log10_NoV), na.rm = TRUE) * 100
  bias  <- mean(val$predicted - val$avg_Log10_NoV, na.rm = TRUE)
  pbias <- (bias / mean(val$avg_Log10_NoV, na.rm = TRUE)) * 100
  corr  <- as.numeric(cor(val$avg_Log10_NoV, val$predicted,
                          use = "complete.obs", method = "spearman"))
  cov   <- COV(val$avg_Log10_NoV,
               lower = val$q0.025,
               upper = val$q0.975)
  
  metrics_fold <- data.frame(
    Fold  = k,   MSE  = mse,  RMSE = rmse,
    MAE   = mae, MAPE = mape, BIAS = bias,
    pBIAS = pbias, CORR = corr, COV = cov
  )
  
  metrics <- rbind(metrics, metrics_fold)
  cat("Fold", k, "complete | RMSE:", round(rmse, 4),
      "| Coverage:", round(cov, 2), "%\n")
}



# ── Overall coverage ──────────────────────────────────────────────────────────
coverage_probability <- mean(
  predictions$avg_Log10_NoV >= predictions$q0.025 &
    predictions$avg_Log10_NoV <= predictions$q0.975,
  na.rm = TRUE
) * 100

cat("Overall Coverage Probability (%):", coverage_probability, "\n")

# ── Summary metrics ───────────────────────────────────────────────────────────
summary_metrics <- metrics %>%
  summarise(
    Mean_MSE       = mean(MSE),
    SE_MSE         = sd(MSE) / sqrt(n()),
    CI_Lower_MSE   = Mean_MSE - qt(0.975, df = n() - 1) * SE_MSE,
    CI_Upper_MSE   = Mean_MSE + qt(0.975, df = n() - 1) * SE_MSE,
    Mean_RMSE      = mean(RMSE),
    SE_RMSE        = sd(RMSE) / sqrt(n()),
    CI_Lower_RMSE  = Mean_RMSE - qt(0.975, df = n() - 1) * SE_RMSE,
    CI_Upper_RMSE  = Mean_RMSE + qt(0.975, df = n() - 1) * SE_RMSE,
    Mean_MAE       = mean(MAE),
    Mean_MAPE      = mean(MAPE),
    Mean_BIAS      = mean(BIAS),
    Mean_pBIAS     = mean(pBIAS),
    Mean_CORR      = mean(CORR),
    COV            = coverage_probability
  )

print(summary_metrics)
