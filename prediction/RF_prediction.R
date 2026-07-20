library(tidyverse)
library(ranger)
library(caret)
library(sf)
library(sp)
library(raster)

set.seed(123)

## read data -------------------------------------------------------------------
nov_df <- read.csv("nov_df_reduced.csv")
pred_grid <- read.csv("lsoa_grid_prediction.csv") 

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
  dplyr::select(avg_Log10_NoV, lockdown_step3, lockdown_step4, lockdown_planB, lockdown_lifting,
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
      scale_school_density = scale_school_den,
      scale_carehome_density = scale_carehome_den,
      scale_mobility = scale_mob_7day,
      scale_BAME = scale_bame,
      scale_imd_score= scale_imd,
      scale_prop_urb= scale_prop_urb,
      scale_rain_rolling_7day = scale_rain_7day_avg,
      scale_temp_rolling_7day = scale_temp_7day_avg,
      date_index = week
    )

grid <- grid %>% 
  dplyr::select(lockdown_step3, lockdown_step4, lockdown_planB, lockdown_lifting,
                scale_school_density, scale_carehome_density, scale_imd_score, scale_BAME,
                scale_mobility, scale_rain_rolling_7day, scale_temp_rolling_7day, scale_prop_urb,Easting, Northing, date_index, LSOA21CD)

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
                    quantreg = TRUE)

# prediction on grid -------------------------------------------------

pred_mean <- predict(model, data = grid)$predictions
pred_quantiles <- predict(model, data = grid, type = "quantiles", quantiles = c(0.025, 0.5, 0.975))$predictions

grid$predicted <- pred_mean
grid$pred_025 <- pred_quantiles[, 1]
grid$pred_50 <- pred_quantiles[, 2]
grid$pred_975 <- pred_quantiles[, 3]

write.csv(grid, "outputs/prediction/rf_lsoa_prediction_new.csv", row.names = FALSE)


