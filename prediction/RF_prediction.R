# =============================================================================
# Norovirus Random Forest Spatio-Temporal Model
# =============================================================================

setwd("~/GitHub/Statistical-Learning-Methods-for-Wastewater-Pathogen-Prediction")

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(ranger)
library(caret)
library(sf)
library(sp)
library(raster)
library(fastDummies)
library(data.table)

set.seed(123)

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

nov_df    <- read.csv("data/nov_df.csv")
pred_grid <- read.csv("data/lsoa_grid_prediction.csv")

# -----------------------------------------------------------------------------
# One-Hot Encode Categorical Variables in Prediction Grid
# -----------------------------------------------------------------------------

pred_grid <- dummy_cols(pred_grid, select_columns = "rain_cat", remove_selected_columns = TRUE)
names(pred_grid) <- sub("^rain_cat_", "rain_", names(pred_grid))

pred_grid <- dummy_cols(pred_grid, select_columns = "temp_cat", remove_selected_columns = TRUE)
names(pred_grid) <- sub("^temp_cat_", "temp_", names(pred_grid))

# -----------------------------------------------------------------------------
# Clean and Scale Covariates
# -----------------------------------------------------------------------------

# Ensure no missing coordinates
coords_master <- nov_df %>%
  filter(!is.na(Easting) & !is.na(Northing)) %>%
  distinct(site_code, Easting, Northing)

nov_df <- nov_df %>%
  dplyr::select(-Easting, -Northing) %>%
  left_join(coords_master, by = "site_code")

covariates_to_scale <- c(
  "school_den", "carehome_den", "imd", "bame",
  "mob_7day", "prop_agri", "prop_urb"
)

scale_covariates <- function(df, covariates_to_scale) {
  scale <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  df %>% mutate(across(all_of(covariates_to_scale), scale, .names = "scale_{.col}"))
}

nov_df <- scale_covariates(nov_df, covariates_to_scale) %>%
  dplyr::select(
    site_code, avg_Log10_NoV,
    lockdown_step3, lockdown_step4, lockdown_planB, lockdown_lifting,
    scale_school_den, scale_carehome_den, scale_imd, scale_bame, scale_mob_7day,
    temp_2, temp_3, temp_4, temp_5,
    rain_2, rain_3, rain_4, rain_5,
    scale_prop_urb, scale_prop_agri,
    Easting, Northing, week
  ) %>%
  filter(!is.na(avg_Log10_NoV)) %>%
  as.data.frame()

# -----------------------------------------------------------------------------
# Model Formula
# -----------------------------------------------------------------------------

formula <- as.formula(
  "avg_Log10_NoV ~ lockdown_step3 + lockdown_step4 + lockdown_planB +
   scale_school_den + scale_carehome_den + scale_mob_7day +
   scale_bame + scale_imd + scale_prop_agri + scale_prop_urb +
   temp_2 + temp_3 + temp_4 + temp_5 +
   rain_2 + rain_3 + rain_4 + rain_5 +
   Easting + Northing + week"
)

# -----------------------------------------------------------------------------
# Train Random Forest Model
# -----------------------------------------------------------------------------

model <- ranger(
  formula,
  data          = nov_df,
  num.trees     = 1000,
  importance    = "permutation",
  mtry          = 8,
  min.node.size = 2,
  write.forest  = TRUE,
  quantreg      = TRUE,
  keep.inbag    = TRUE
)

# -----------------------------------------------------------------------------
# Predict on LSOA Grid
# -----------------------------------------------------------------------------

predictions   <- predict(model, data = pred_grid, predict.all = TRUE)
pred_quantiles <- predict(model, data = pred_grid, type = "quantiles", quantiles = c(0.025, 0.5, 0.975))$predictions
pred_sd        <- predict(model, data = pred_grid, type = "se", se.method = "infjack", keep.inbag = TRUE)$se

pred_grid <- pred_grid %>%
  dplyr::select(LSOA21CD, week) %>%
  mutate(
    mean  = rowMeans(predictions$predictions),
    lower = pred_quantiles[, 1],
    upper = pred_quantiles[, 3],
    sd    = pred_sd
  )

# -----------------------------------------------------------------------------
# Save Mean Predictions
# -----------------------------------------------------------------------------

write.csv(pred_grid, "outputs/rf_lsoa_prediction_mean.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Attach Posterior Samples and Save
# -----------------------------------------------------------------------------

pred_samples <- as.data.frame(predictions$predictions)
colnames(pred_samples) <- paste0("sample_", 1:ncol(pred_samples))

pred_grid <- cbind(pred_grid, pred_samples)

fwrite(pred_grid, "outputs/rf_lsoa_prediction_samples.csv")

# -----------------------------------------------------------------------------
# National Weekly Predictions (Population-Weighted)
# -----------------------------------------------------------------------------

lsoa_pop <- read.csv("data/population_lsoa_21.csv")

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

write.csv(nat_pred, "outputs/rf_nat_pred.csv", row.names = FALSE)

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

lsoa_sd <- pred_grid %>%
  group_by(LSOA21CD) %>%
  summarise(sd = mean(sd), .groups = "drop")

sample_mat <- as.matrix(lsoa_samples[, ..sample_cols])

lsoa_summary <- data.frame(
  LSOA21CD = lsoa_samples$LSOA21CD,
  mean     = rowMeans(sample_mat),
  lower    = apply(sample_mat, 1, quantile, 0.025),
  upper    = apply(sample_mat, 1, quantile, 0.975),
  sd       = lsoa_sd$sd
)

write.csv(lsoa_summary, "outputs/rf_lsoa_pred.csv", row.names = FALSE)