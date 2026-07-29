# =============================================================================
# Norovirus INLA Spatio-Temporal Model
# =============================================================================

setwd("~/GitHub/Statistical-Learning-Methods-for-Wastewater-Pathogen-Prediction")

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(INLA)
library(inlabru)
library(sp)
library(fmesher)
library(corrr)
library(sf)
library(ggplot2)
library(viridis)
library(ggforce)
library(reshape2)
library(readr)
library(fastDummies)
library(data.table)

inla.setOption(num.threads = "6")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

nov_df   <- read.csv("data/nov_df.csv")
pred_grid <- fread("data/lsoa_grid_prediction.csv")

england <- st_read("data/england_crop.shp")
england <- st_transform(england, crs = 27700)

# -----------------------------------------------------------------------------
# Scale Covariates
# -----------------------------------------------------------------------------

covariates_to_scale <- c(
  "school_den", "carehome_den", "imd", "bame",
  "mob_7day", "prop_agri", "prop_urb"
)

scale_covariates <- function(df, covariates_to_scale) {
  scale <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  df %>% mutate(across(all_of(covariates_to_scale), scale, .names = "scale_{.col}"))
}

nov_df <- scale_covariates(nov_df, covariates_to_scale) %>% as.data.frame()

# Rebuild coordinates from master to avoid duplicates/NAs
coords_master <- nov_df %>%
  filter(!is.na(Easting) & !is.na(Northing)) %>%
  distinct(site_code, Easting, Northing)

nov_df <- nov_df %>%
  dplyr::select(-Easting, -Northing) %>%
  left_join(coords_master, by = "site_code")

# Convert to sf and reproject to km units
bng_km <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

nov <- st_as_sf(nov_df, coords = c("Easting", "Northing"), crs = 27700)
nov <- st_transform(nov, crs = bng_km)

# -----------------------------------------------------------------------------
# One-Hot Encode Categorical Variables in Prediction Grid
# -----------------------------------------------------------------------------

pred_grid <- dummy_cols(pred_grid, select_columns = "rain_cat", remove_selected_columns = TRUE)
names(pred_grid) <- sub("^rain_cat_", "rain_", names(pred_grid))

pred_grid <- dummy_cols(pred_grid, select_columns = "temp_cat", remove_selected_columns = TRUE)
names(pred_grid) <- sub("^temp_cat_", "temp_", names(pred_grid))

# -----------------------------------------------------------------------------
# Priors and Model Formula
# -----------------------------------------------------------------------------

pc_prec         <- list(prec = list(prior = "pc.prec", param = c(10, 0.05)))
rho_hyper       <- list(theta = list(prior = "pccor0", param = c(0.1, 0.9)))

formula <- as.formula(
  'avg_Log10_NoV ~ -1 + Intercept +
   lockdown_step3 + lockdown_step4 + lockdown_planB +
   scale_school_den + scale_carehome_den + scale_mob_7day +
   scale_bame + scale_imd + scale_prop_agri + scale_prop_urb +
   temp_2 + temp_3 + temp_4 + temp_5 +
   rain_2 + rain_3 + rain_4 + rain_5 +
   f(site_code, model = "iid", hyper = pc_prec) +
   f(week, model = "iid", hyper = pc_prec) +
   f(spatial.field, model = spde, group = spatial.field.group,
     control.group = list(model = "iid", hyper = pc_prec))'
)

# -----------------------------------------------------------------------------
# Build Mesh
# -----------------------------------------------------------------------------

coords       <- as.matrix(st_coordinates(nov))
max.edge     <- 10
bound.outer  <- diff(range(coords[, 1])) / 5
domain       <- inla.nonconvex.hull(coords)

mesh <- fm_mesh_2d_inla(
  boundary = domain,
  loc      = coords,
  max.edge = c(1, 2) * max.edge,
  offset   = c(max.edge, bound.outer),
  cutoff   = max.edge / 5,
  crs      = st_crs(nov)
)

# -----------------------------------------------------------------------------
# SPDE Object
# -----------------------------------------------------------------------------

spde <- inla.spde2.pcmatern(
  mesh        = mesh,
  alpha       = 2,
  prior.range = c(10, 0.05),
  prior.sigma = c(5, 0.05)
)

# -----------------------------------------------------------------------------
# Estimation Stack
# -----------------------------------------------------------------------------

n_week  <- length(unique(nov$week))
s.index <- inla.spde.make.index(
  name    = "spatial.field",
  group   = nov$week,
  n.spde  = spde$n.spde,
  n.group = n_week
)

A_est <- inla.spde.make.A(
  mesh    = mesh,
  loc     = coords,
  group   = nov$week,
  n.group = n_week
)

stack_est <- inla.stack(
  tag  = "est",
  data = list(avg_Log10_NoV = nov$avg_Log10_NoV),
  A    = list(A_est, 1),
  effects = list(
    c(s.index, list(Intercept = 1)),
    list(
      week              = nov$week,
      site_code         = nov$site_code,
      lockdown_step3    = nov$lockdown_step3,
      lockdown_step4    = nov$lockdown_step4,
      lockdown_lifting  = nov$lockdown_lifting,
      lockdown_planB    = nov$lockdown_planB,
      scale_school_den  = nov$scale_school_den,
      scale_carehome_den = nov$scale_carehome_den,
      scale_mob_7day    = nov$scale_mob_7day,
      scale_bame        = nov$scale_bame,
      scale_imd         = nov$scale_imd,
      scale_prop_urb    = nov$scale_prop_urb,
      scale_prop_agri   = nov$scale_prop_agri,
      temp_2 = nov$temp_2, temp_3 = nov$temp_3,
      temp_4 = nov$temp_4, temp_5 = nov$temp_5,
      rain_2 = nov$rain_2, rain_3 = nov$rain_3,
      rain_4 = nov$rain_4, rain_5 = nov$rain_5
    )
  )
)

# -----------------------------------------------------------------------------
# Prediction Stack
# -----------------------------------------------------------------------------

sc           <- 1 / 1000
coords.lsoa  <- as.matrix(cbind(pred_grid[, "Easting"], pred_grid[, "Northing"])) * sc
n_week_pred  <- length(unique(pred_grid$week))
group_pred   <- pred_grid$week

A_pred <- inla.spde.make.A(
  mesh    = mesh,
  loc     = coords.lsoa,
  group   = group_pred,
  n.group = n_week_pred
)

s.index_pred <- inla.spde.make.index(
  name    = "spatial.field",
  n.spde  = spde$n.spde,
  n.group = n_week_pred
)

stack_pred <- inla.stack(
  tag  = "pred",
  data = list(avg_Log10_NoV = NA),
  A    = list(A_pred, 1),
  effects = list(
    c(s.index_pred, list(Intercept = 1)),
    list(
      week               = pred_grid$week,
      LSOA21CD           = pred_grid$LSOA21CD,
      lockdown_step3     = pred_grid$lockdown_step3,
      lockdown_step4     = pred_grid$lockdown_step4,
      lockdown_lifting   = pred_grid$lockdown_lifting,
      lockdown_planB     = pred_grid$lockdown_planB,
      scale_school_den   = pred_grid$scale_school_den,
      scale_carehome_den = pred_grid$scale_carehome_den,
      scale_mob_7day     = pred_grid$scale_mob_7day,
      scale_bame         = pred_grid$scale_bame,
      scale_imd          = pred_grid$scale_imd,
      scale_prop_urb     = pred_grid$scale_prop_urb,
      scale_prop_agri    = pred_grid$scale_prop_agri,
      temp_2 = pred_grid$temp_2, temp_3 = pred_grid$temp_3,
      temp_4 = pred_grid$temp_4, temp_5 = pred_grid$temp_5,
      rain_2 = pred_grid$rain_2, rain_3 = pred_grid$rain_3,
      rain_4 = pred_grid$rain_4, rain_5 = pred_grid$rain_5
    )
  )
)

# -----------------------------------------------------------------------------
# Join Stacks
# -----------------------------------------------------------------------------

join_stack <- inla.stack(stack_est, stack_pred)

# -----------------------------------------------------------------------------
# Fit INLA Model
# -----------------------------------------------------------------------------

pred.lsoa <- inla(
  formula,
  data             = inla.stack.data(join_stack),
  family           = "gaussian",
  control.predictor = list(
    A       = inla.stack.A(join_stack),
    compute = TRUE,
    link    = 1
  ),
  control.compute = list(
    dic                          = TRUE,
    cpo                          = TRUE,
    waic                         = TRUE,
    mlik                         = TRUE,
    return.marginals.predictor   = TRUE,
    config                       = TRUE,
    openmp.strategy              = "default",
    smtp                         = "taucs"
  ),
  control.fixed = list(
    mean            = 0,
    prec            = 0.0001,
    mean.intercept  = 0,
    prec.intercept  = 0.0001
  ),
  control.family = list(hyper = list(prec = list(param = c(1, 0.2)))),
  verbose = FALSE
)

# -----------------------------------------------------------------------------
# Extract Predictions and Save
# -----------------------------------------------------------------------------

index <- inla.stack.index(join_stack, "pred")$data

mean_pred <- data.frame(
  LSOA21CD  = pred_grid$LSOA21CD,
  week      = pred_grid$week,
  pred_mean = pred.lsoa$summary.fitted.values[index, "mean"],
  pred_ll   = pred.lsoa$summary.fitted.values[index, "0.025quant"],
  pred_ul   = pred.lsoa$summary.fitted.values[index, "0.975quant"],
  sd        = pred.lsoa$summary.fitted.values[index, "sd"]
)

fwrite(mean_pred, "outputs/inla_prediction.csv")

# -----------------------------------------------------------------------------
# Draw Posterior Samples and Save
# -----------------------------------------------------------------------------

posterior.samples <- inla.posterior.sample(
  n                 = 1000,
  result            = pred.lsoa,
  intern            = FALSE,
  use.improved.mean = TRUE,
  add.names         = TRUE,
  seed              = 0L
)

# Extract fitted values for prediction locations across all samples
samples.fitted <- matrix(
  unlist(lapply(
    lapply(posterior.samples, "[[", "latent"),
    "[", inla.stack.index(join_stack, "pred")$data
  ), use.names = FALSE),
  nrow    = nrow(pred_grid),
  byrow   = FALSE
)

samples.fitted <- as.data.frame(samples.fitted)
colnames(samples.fitted)[1:100] <- paste0("sample_", 1:1000)

samples.fitted <- samples.fitted %>%
  mutate(
    LSOA21CD = pred_grid$LSOA21CD,
    week     = pred_grid$week,
    mean     = pred.lsoa$summary.fitted.values[index, "mean"],
    lower    = pred.lsoa$summary.fitted.values[index, "0.025quant"],
    upper    = pred.lsoa$summary.fitted.values[index, "0.975quant"],
    sd       = pred.lsoa$summary.fitted.values[index, "sd"]
  )

# Convert sample columns to numeric and save
sample_cols <- paste0("sample_", 1:1000)
samples.fitted[, sample_cols] <- apply(samples.fitted[, sample_cols], 2, as.numeric)

fwrite(samples.fitted, "outputs/inla_lsoa_prediction_samples.csv")

# -----------------------------------------------------------------------------
# National Weekly Predictions (Population-Weighted)
# -----------------------------------------------------------------------------

pred_grid <- fread("outputs/inla_lsoa_prediction_samples.csv")
lsoa_pop  <- read.csv("data/population_lsoa_21.csv")

setDT(pred_grid)
setDT(lsoa_pop)

sample_cols <- grep("^sample_", names(pred_grid), value = TRUE)

weekly_samples <- lsoa_pop[pred_grid, on = "LSOA21CD"][
  ,
  lapply(.SD, function(x) weighted.mean(x, w = population, na.rm = TRUE)),
  by      = week,
  .SDcols = sample_cols
]

sample_mat <- as.matrix(weekly_samples[, ..sample_cols])

nat_pred <- data.frame(
  week  = weekly_samples$week,
  mean  = rowMeans(sample_mat),
  lower = apply(sample_mat, 1, quantile, 0.025),
  upper = apply(sample_mat, 1, quantile, 0.975)
)

write.csv(nat_pred, "outputs/inla_nat_pred.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# LSOA-Level Summary (Mean Across Weeks)
# -----------------------------------------------------------------------------

sample_cols <- grep("^sample_", names(pred_grid), value = TRUE)

lsoa_samples <- pred_grid[
  ,
  lapply(.SD, mean, na.rm = TRUE),
  by      = LSOA21CD,
  .SDcols = sample_cols
]

sample_mat <- as.matrix(lsoa_samples[, ..sample_cols])

lsoa_sd<- mean_pred %>%
  group_by(LSOA21CD) %>%
  summarise(sd= mean(sd))

lsoa_summary <- data.frame(
  LSOA21CD = lsoa_samples$LSOA21CD,
  mean     = rowMeans(sample_mat),
  lower    = apply(sample_mat, 1, quantile, 0.025),
  upper    = apply(sample_mat, 1, quantile, 0.975),
  sd= lsoa_sd$sd
)

write.csv(lsoa_summary, "outputs/inla_lsoa_pred.csv", row.names = FALSE)