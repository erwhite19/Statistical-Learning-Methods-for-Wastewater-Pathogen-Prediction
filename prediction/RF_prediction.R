library(tidyverse)
library(ranger)
library(caret)
library(sf)
library(sp)
library(raster)
library(fastDummies)
set.seed(123)


## read data -------------------------------------------------------------------
nov_df <- read.csv("data/nov_df.csv")
pred_grid <- read.csv("data/lsoa_grid_prediction.csv") 


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
  
  
pred_grid$mean <- pred_mean
pred_grid$lower <- pred_quantiles[, 1]
pred_grid$upper <- pred_quantiles[, 3]
pred_grid$sd <- pred_sd

pred_grid<- pred_grid %>%
  dplyr::select(LSOA21CD, week, mean, lower, upper, sd)


# save
write.csv(pred_grid, "outputs/rf_lsoa_prediction_mean.csv", row.names = FALSE)

pred_samples<- predictions$predictions # Rename columns of the predictions matrix

colnames(pred_samples) <- paste0("sample_", 1:ncol(predictions$predictions))

pred_grid<- cbind(pred_grid, pred_samples)

fwrite(pred_grid, "output/rf_lsoa_prediction_samples.csv", row.names = FALSE)

sample_cols <- grep("^sample_", names(pred_grid), value = TRUE)

# lsoa population
lsoa_pop<- read.csv("data/population_lsoa_21.csv")


# Get weekly mean (use data.table to help with large dataset)

# Convert both to data.table (does this in place, no copy made)
setDT(pred_grid)
setDT(lsoa_pop)




# Join + weighted mean in data.table
weekly_samples <- lsoa_pop[pred_grid, on = "LSOA21CD"] [  # left join
  ,                                                         # no row filter
  lapply(.SD, function(x) weighted.mean(x, w = population, na.rm = TRUE)),
  by = week,                                                # group by week
  .SDcols = sample_cols                                     # only sample columns
]


# compute intervals
sample_cols <- grep("^sample_", names(weekly_samples), value = TRUE)
sample_mat <- as.matrix(weekly_samples[, ..sample_cols])

nat_pred<- data.frame(
  week     = weekly_samples$week,
  mean = rowMeans(sample_mat),
  lower  = apply(sample_mat, 1, quantile, 0.025),
  upper  = apply(sample_mat, 1, quantile, 0.975)
)


write.csv(weekly_summary, "output/rf_nat_pred.csv")

# Get LSOA mean (use data.table to help with large dataset)

lsoa_samples <- pred_grid[  # left join
  ,                                                         # no row filter
  lapply(.SD, function(x) mean(x, w = population, na.rm = TRUE)),
  by = LSOA21CD,                                                # group by week
  .SDcols = sample_cols                                     # only sample columns
]

lsoa_sd<- pred_grid %>%
  group_by(LSOA21CD) %>%
  summarise(sd= mean(sd))

# compute intervals
sample_cols <- grep("^sample_", names(lsoa_samples), value = TRUE)
sample_mat <- as.matrix(lsoa_samples[, ..sample_cols])

lsoa_summary <- data.frame(
  LSOA21CD     = lsoa_samples$LSOA21CD,
  mean = rowMeans(sample_mat),
  lower  = apply(sample_mat, 1, quantile, 0.025),
  upper  = apply(sample_mat, 1, quantile, 0.975),
  sd= lsoa_sd$sd
)


write.csv(lsoa_summary, "outputs/rf_lsoa_pred.csv")
