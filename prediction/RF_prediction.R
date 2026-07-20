library(tidyverse)
library(ranger)
library(caret)
library(sf)
library(sp)
library(raster)
library(fastDummies)
set.seed(123)
setwd("C:/Users/erw19/OneDrive - Imperial College London/PhD/script")

## read data -------------------------------------------------------------------
nov_df <- read.csv("nov_df_reduced.csv")
pred_grid <- read.csv("Data/prediction_data/lsoa_grid_prediction.csv") 


pred_grid <- dummy_cols(
  pred_grid,
  select_columns = "rain_cat",
  remove_selected_columns = TRUE
)

names(pred_grid) <- sub("^rain_cat_", "rain_", names(pred_grid))

pred_grid <- dummy_cols(
  pred_grid,
  select_columns = "temp_cat",
  remove_selected_columns = TRUE
)

names(pred_grid) <- sub("^temp_cat_", "temp_", names(pred_grid))


# ensure no missing coordinates
coords_master <- nov_df %>%
  filter(!is.na(Easting) & !is.na(Northing)) %>%
  distinct(site_code, Easting, Northing)

nov_df <- nov_df %>%
  dplyr::select(-Easting, -Northing) %>%
  left_join(coords_master, by = "site_code")

nov_df <- nov_df %>%
  mutate(avg_Log10_NoV = ifelse(GI_ND | GII_ND, NA, avg_Log10_NoV)) 

## scale------------------------------------------------------------------------

covariates_to_scale <-  c("school_den", "carehome_den", "imd", "bame", "mob_7day", "prop_agri", "prop_urb")
scale_covariates <- function(df, covariates_to_scale) {
  # Define a custom scaling function
  scale <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  # Apply the custom scaling function across specified covariates
  df <- df %>% mutate(across(all_of(covariates_to_scale), scale, .names = "scale_{.col}"))
  
  return(df)
}

nov_df<- scale_covariates(nov_df, covariates_to_scale)

nov_df <- nov_df %>% 
  dplyr::select(site_code, avg_Log10_NoV, lockdown_step3, lockdown_step4, lockdown_planB, lockdown_lifting,
                scale_school_den, scale_carehome_den, scale_imd, scale_bame,
                scale_mob_7day,temp_2, temp_3, temp_4, temp_5, rain_2, rain_3 , rain_4 , rain_5, scale_prop_urb,scale_prop_agri,
                Easting, Northing, week)


nov_df<- as.data.frame(nov_df)

#remove rows with NA in target variable
nov_df <- nov_df[!is.na(nov_df$avg_Log10_NoV), ]


# convert df into sf and units m to km
nov <- st_as_sf(nov_df, coords= c("Easting", "Northing"), crs= 27700)
nov <- st_transform(nov,  crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")

#grid <-st_as_sf(pred_grid, coords= c("Easting", "Northing"), crs= 27700)
#grid <- st_transform(grid,  crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")

# rename columns in pred_grid -------------------------------------------------

grid <- pred_grid %>%
  rename(lockdown_step3 = lockdown_step3,
      lockdown_step4 = lockdown_step4,
      lockdown_lifting = lockdown_lifting,
      lockdown_planB = lockdown_planB,
      scale_school_den = scale_school_den,
      scale_carehome_density = scale_carehome_den,
      scale_mobility = scale_mob_7day,
      scale_ = scale_bame,
      scale_imd_score= scale_imd,
      scale_prop_urb= scale_prop_urb,
      scale_rain_rolling_7day = scale_rain_7day_avg,
      scale_temp_rolling_7day = scale_temp_7day_avg,
      week = week
    )


# formula ---------------------------------------------------------------
formula <- as.formula("avg_Log10_NoV ~ lockdown_step3 + lockdown_step4 + lockdown_planB +
    scale_school_den + scale_carehome_den + scale_mob_7day + scale_bame + scale_imd + scale_prop_agri + scale_prop_urb +
    temp_2 + temp_3 + temp_4 + temp_5 +
    rain_2 + rain_3 + rain_4 + rain_5+ Easting + Northing + week")

# train model----------------------------------------------------------------

model  <- ranger(formula,
                    data=nov_df, 
                    num.trees = 1000, 
                    importance = "permutation",
                    mtry = 8, 
                    min.node.size=2,
                    write.forest = TRUE,
                    quantreg = TRUE,
                 keep.inbag = T)

# prediction on grid -------------------------------------------------

predictions <- predict(model, data = pred_grid, predict.all= T)

pred_mean <- rowMeans(predictions$predictions)

# Uncertainty

pred_quantiles <- predict(model, data = pred_grid, type = "quantiles", quantiles = c(0.025, 0.5, 0.975))$predictions

pred_sd<- predict(model, data = pred_grid, type = "se", se.method = "infjack", keep.inbag= T)$se
  
  
pred_grid$predicted <- pred_mean
pred_grid$pred_025 <- pred_quantiles[, 1]
pred_grid$pred_50 <- pred_quantiles[, 2]
pred_grid$pred_975 <- pred_quantiles[, 3]
pred_grid$sd <- pred_sd

# save
write.csv(pred_grid, ".../rf_lsoa_prediction.csv", row.names = FALSE)

pred_samples<- predictions$predictions # Rename columns of the predictions matrix

colnames(pred_samples) <- paste0("sample_", 1:ncol(predictions$predictions))

pred_grid<- cbind(pred_grid, pred_samples)

sample_cols <- grep("^sample_", names(pred_grid), value = TRUE)

# Get weekly mean for EACH sample separately
weekly_samples <- pred_grid %>%
  group_by(week) %>%
  summarise(across(all_of(sample_cols), mean), .groups = "drop")

# compute intervals
sample_mat <- as.matrix(weekly_samples[, sample_cols])

weekly_summary <- data.frame(
  week     = weekly_samples$week,
  mean = rowMeans(sample_mat),
  lower  = apply(sample_mat, 1, quantile, 0.025),
  upper  = apply(sample_mat, 1, quantile, 0.975)
)


write.csv(weekly_summary, ".../rf_nat_pred.csv")


