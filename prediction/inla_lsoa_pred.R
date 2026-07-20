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

#setwd("~/Term3-project")

inla.setOption(num.threads= "6")

## read data -------------------------------------------------------------------
nov_df <- read.csv("nov_df_reduced.csv")
pred_grid <- read.csv("Data/prediction_data/lsoa_grid_prediction.csv")

#englang shapefile
england<- st_read("Data/shapefiles/england/england_crop.shp")
england<- st_transform(england, crs= 27700)

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
nov_df<- as.data.frame(nov_df)

coords_master <- nov_df %>%
  filter(!is.na(Easting) & !is.na(Northing)) %>%
  distinct(site_code, Easting, Northing)

nov_df <- nov_df %>%
  dplyr::select(-Easting, -Northing) %>%
  left_join(coords_master, by = "site_code")

# convert df into sf and units m to km
nov <- st_as_sf(nov_df, coords= c("Easting", "Northing"), crs= 27700)
nov <- st_transform(nov,  crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")

#grid <-st_as_sf(pred_grid, coords= c("Easting", "Northing"), crs= 27700)
#grid <- st_transform(grid,  crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")

#one hot encode

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



## Set prior and formula--------------------------------------------------------

prec_overall_time <- c(10, 0.05)
prec_overall_space <- c(10, 0.05)

intercept_prior <- list(mean.intercept = 0, prec.intercept = 0.0001)
fixed_effects_prior <- list(mean = 0, prec = 0.0001)


rho_hyper = list(theta = list(prior = "pccor0", param = c(0.1, 0.9)))
pc_prec <- list(prec = list(prior = "pc.prec", param = c(10, 0.05)))

formula<-  as.formula('avg_Log10_NoV ~ -1 + Intercept + lockdown_step3 + lockdown_step4 + lockdown_planB +
    scale_school_den + scale_carehome_den + scale_mob_7day + scale_bame + scale_imd + scale_prop_agri + scale_prop_urb +
    temp_2 + temp_3 + temp_4 + temp_5 +
    rain_2 + rain_3 + rain_4 + rain_5 +
                        f(site_code, model="iid", hyper= pc_prec) + f(week, model= "iid", hyper= pc_prec) +
                        f(spatial.field, model=spde, group=spatial.field.group, control.group=list(model="iid", hyper=pc_prec))')



## mesh(est) ------------------------------------------------------------------------

max.edge = 10
bound.outer = diff(range(st_coordinates(nov)[,1]))/5
coords<- as.matrix((st_coordinates(nov)))
domain <- inla.nonconvex.hull(coords)

mesh<- fm_mesh_2d_inla(boundary = domain,
                       loc=coords,
                       max.edge = c(1,2)*max.edge,
                       offset= c(max.edge, bound.outer),
                       cutoff= max.edge/5,
                       crs= st_crs(nov))

## SPDE(est) -------------------------------------------------------------------------

spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2,
                            prior.range = c(10, 0.05),
                            prior.sigma = c(5, 0.05))

# create index
n_week<- length(unique(nov$week))
s.index<- inla.spde.make.index(name="spatial.field",
                               group = nov$week,
                               n.spde=spde$n.spde,
                               n.group= n_week) 
## A matrix (est)---------------------------------------------------------------

A_est <- inla.spde.make.A(mesh = mesh,
                          loc = coords,
                          group = nov$week,
                          n.group = n_week)

dim(A_est)

## estimation stack ------------------------------------------------------------

stack_est <- inla.stack(
  tag = "est",
  data = list(avg_Log10_NoV = nov$avg_Log10_NoV),
  A = list(A_est,1),
  effects = list(
    c(s.index, list(Intercept = 1)),
    list(
      week = nov$week,
      site_code = nov$site_code,
      lockdown_step3 = nov$lockdown_step3,
      lockdown_step4 = nov$lockdown_step4,
      lockdown_lifting = nov$lockdown_lifting,
      lockdown_planB = nov$lockdown_planB,
      scale_school_den = nov$scale_school_den,
      scale_carehome_den = nov$scale_carehome_den,
      scale_mob_7day = nov$scale_mob_7day,
      scale_bame = nov$scale_bame,
      scale_imd= nov$scale_imd,
      scale_prop_urb= nov$scale_prop_urb,
      scale_prop_agri = nov$scale_prop_agri,
      temp_2 = nov$temp_2,
      temp_3 = nov$temp_3,
      temp_4 = nov$temp_4,
      temp_5 = nov$temp_5,
      rain_2 = nov$rain_2,
      rain_3 = nov$rain_3,
      rain_4 = nov$rain_4,
      rain_5 = nov$rain_5
    )
  )
)

## prediction stack ------------------------------------------------------------

sc <- 1/1000
coords.lsoa<- as.matrix(cbind(pred_grid[,"Easting"], pred_grid[,"Northing"])*sc)
#coords.lsoa<- st_coordinates(grid)

nrow(coords.lsoa) == nrow(grid)


group_pred <- pred_grid$week
n_week_pred <- length(unique(pred_grid$week))

## A matrix (pred) -------------------------------------------------------------

A_pred = inla.spde.make.A(
  mesh =  mesh,
  loc = coords.lsoa,
  group = group_pred,
  n.group = n_week_pred)

dim(A_pred)

## create index (pred) ---------------------------------------------------------

s.index_pred = inla.spde.make.index(
  name = "spatial.field",
  n.spde = spde$n.spde,
  n.group = n_week_pred)

## prediction stack-------------------------------------------------------------

stack_pred <- inla.stack(
  tag = "pred",
  data = list(avg_Log10_NoV = NA),
  A = list(A_pred,1),
  effects = list(
    c(s.index_pred, list(Intercept = 1)),
    list(
      week = pred_grid$week,
      LSOA21CD = pred_grid$LSOA21CD,
      lockdown_step3 = pred_grid$lockdown_step3,
      lockdown_step4 = pred_grid$lockdown_step4,
      lockdown_lifting = pred_grid$lockdown_lifting,
      lockdown_planB = pred_grid$lockdown_planB,
      scale_school_den = pred_grid$scale_school_den,
      scale_carehome_den = pred_grid$scale_carehome_den,
      scale_mob_7day = pred_grid$scale_mob_7day,
      scale_bame = pred_grid$scale_bame,
      scale_imd= pred_grid$scale_imd,
      scale_prop_urb= pred_grid$scale_prop_urb,
      scale_prop_agri = pred_grid$scale_prop_agri,
      temp_2 = pred_grid$temp_2,
      temp_3 = pred_grid$temp_3, 
      temp_4 = pred_grid$temp_4,
      temp_5 = pred_grid$temp_5,
      rain_2 = pred_grid$rain_2,
      rain_3 = pred_grid$rain_3,
      rain_4 = pred_grid$rain_4,
      rain_5 = pred_grid$rain_5
     
    )
  )
)

## join stack ------------------------------------------------------------------

join_stack <- inla.stack(stack_est, stack_pred)

## prediction model 

pred.lsoa <- inla(formula,
                  data = inla.stack.data(join_stack),
                  family = "gaussian",
                  control.predictor = list(A = inla.stack.A(join_stack), compute = TRUE, link = 1),
                  control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE,
                                         mlik = TRUE, return.marginals.predictor = TRUE, config = TRUE,
                                         openmp.strategy = "default", smtp = "taucs"),
                  control.fixed = list(mean = 0, prec = 0.0001,mean.intercept = 0, prec.intercept = 0.0001), 
                  control.family = list(hyper = list(prec = list(param = c(1, 0.2)))),
                  verbose= F)



index <- inla.stack.index(join_stack, "pred")$data

# mean predeictions

mean_pred <- data.frame(
  LSOA21CD = pred_grid$LSOA21CD,
  week      = pred_grid$week,
  pred_mean = pred.lsoa$summary.fitted.values[index, "mean"],
  pred_ll   = pred.lsoa$summary.fitted.values[index, "0.025quant"],
  pred_ul   = pred.lsoa$summary.fitted.values[index, "0.975quant"],
  sd        = pred.lsoa$summary.fitted.values[index, "sd"]
)

fwrite(mean_pred, "output/inla_mean_predictions.csv")
# extract sample.size values from each marginal of the fitted values
posterior.samples = inla.posterior.sample(n = 100, result=pred.lsoa, intern = FALSE, use.improved.mean = TRUE, add.names = TRUE, seed = 0L)

# extract the posterior latent for the validation sites and reshape in order to have a matrix of nrow(valid) x sample.size
samples.fitted<- matrix(
  unlist(lapply(
    lapply(posterior.samples, "[[", "latent"), 
    "[", inla.stack.index(join_stack, "pred")$data
  ), use.names = FALSE), 
  nrow = nrow(pred_grid), 
  byrow = FALSE
)

samples.fitted <- as.data.frame(samples.fitted) %>%
  mutate(
    LSOA21CD = pred_grid$LSOA21CD,
    week      = pred_grid$week,
    pred_ll   = pred.lsoa$summary.fitted.values[index, "0.025quant"],
    pred_ul   = pred.lsoa$summary.fitted.values[index, "0.975quant"],
    sd        = pred.lsoa$summary.fitted.values[index, "sd"]
  )


# Create the column names for columns 1 to 100
new_names <- paste0("sample_", 1:100)

# Assign the new names to columns 1 to 100
colnames(samples.fitted)[1:100] <- new_names
columns_to_convert <- 1:100

# Step 3: Convert specified columns to numeric
samples.df<- as.data.frame(samples.fitted)
samples.df[, columns_to_convert] <- apply(samples.df[, columns_to_convert], 2, as.numeric)

fwrite(samples.df, "output/inla_predictions.csv")



## Aggregate by time for plotting

time_agg <- samples.df %>%
  group_by(week) %>%
  summarise(across(as.numeric(starts_with("sample_")), \(x) mean(x, na.rm = TRUE)))




# Function to calculate mean and 95% confidence interval
calculate_stats <- function(x) {
  mean_val <- mean(x)
  lower <- quantile(x, 0.025)
  upper <- quantile(x, 0.975)
  data.frame(mean = mean_val, lower = lower, upper = upper)
}

# Gather the data into long format, perform calculations, and spread it back
time_long <- time_agg %>%
  pivot_longer(cols = starts_with("sample_"), names_to = "sample", values_to = "value") %>%
  group_by(week) %>%
  summarise(
    calculate_stats(value),
    .groups = 'drop'
  )

# save outputs
write.csv(time_long, "output/inla_national_predictions.csv") # aggregated predictions

space_agg<- samples.df %>%
  group_by(LSOA21CD) %>%
  summarise(across(as.numeric(starts_with("sample_")), \(x) mean(x, na.rm = TRUE)))

# Gather the data into long format, perform calculations, and spread it back
space_long <- space_agg %>%
  pivot_longer(cols = starts_with("sample_"), names_to = "sample", values_to = "value") %>%
  group_by(LSOA21CD) %>%
  summarise(
    calculate_stats(value),
    .groups = 'drop'
  )



write.csv(space_long, "output/inla_lsoa_mean.csv") # aggregated predictions


if(FALSE){


samples.df<- fread("outputs/wales/NoV/predictions/m1_predictions.csv")


results <- samples.df %>%
  pivot_longer(cols = starts_with("sample_"), names_to = "draw", values_to = "value") %>%
  group_by(LSOA21CD, week) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  )

# replicate so we have n lsoas x 46 
study_start <- as.Date("2022-09-05")

weeks <- 1:46  # Number of study weeks

results<- results %>%
  mutate(week_start = study_start + (week - 1) * 7) %>% 
  arrange(LSOA21CD, week)

write.csv(results, "outputs/wales/NoV/predictions/lsoa_week_preds.csv")

#overall average for LSOA prediction
overall_avg <- pred_grid %>%
  group_by(LSOA21CD, Easting, Northing) %>%
  summarise(pred_mean = mean(pred_mean), .groups = "drop")

overall_sd <- pred_grid %>%
  group_by(LSOA21CD, Easting, Northing) %>%
  summarise(sd = mean(sd), .groups = "drop")


# LSOA boundaries
lsoa_boundaries <- st_read("Data/LSOA/LSOA2021_boundaries/LSOA2021_boundaries.shp") %>%
  dplyr::select(LSOA21CD, geometry)

lsoa_boundaries<- lsoa_boundaries[grepl("^E0", lsoa_boundaries$LSOA21CD), ]

lsoa_joined_pred <- lsoa_boundaries %>%
  left_join(overall_avg, by = c("LSOA21CD" = "LSOA21CD")) %>%
  left_join(overall_sd, by = c("LSOA21CD" = "LSOA21CD")) %>%
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700)

pdf("lsoa_prediction_plot.pdf",width = 10, height =8)
ggplot(lsoa_joined_pred) +
  geom_sf(aes(fill = pred_mean), color = NA) +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "LSOA prediction: INLA",
       fill = "NoV concentration (Log10gc per 100k capita)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

  dev.off()
  
  pdf("lsoa_prediction_plot_sd.pdf",width = 10, height =8)
  ggplot(lsoa_joined_pred) +
    geom_sf(aes(fill = sd), color = NA) +
    scale_fill_viridis_c(option = "mako") +
    labs(title = "LSOA prediction: INLA",
         fill = "Standard deviation") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  dev.off()
  
pred_grid_rf<- read.csv("output/rf_lsoa_prediction_new.csv")

#overall average for LSOA prediction
overall_avg <- pred_grid_rf %>%
  group_by(LSOA21CD, Easting, Northing) %>%
  summarise(pred_mean = mean(predicted), .groups = "drop")

overall_sd <- pred_grid_rf %>%
  group_by(LSOA21CD, Easting, Northing) %>%
  summarise(sd = mean(sd), .groups = "drop")

  
  # LSOA boundaries
  lsoa_boundaries <- st_read("Data/LSOA/LSOA2021_boundaries/LSOA2021_boundaries.shp") %>%
    dplyr::select(LSOA21CD, geometry)
  
  lsoa_boundaries<- lsoa_boundaries[grepl("^E0", lsoa_boundaries$LSOA21CD), ]
  
  lsoa_joined_pred <- lsoa_boundaries %>%
    left_join(overall_avg, by = c("LSOA21CD" = "LSOA21CD")) %>%
    left_join(overall_sd, by = c("LSOA21CD" = "LSOA21CD")) %>%
    st_as_sf(coords = c("Easting", "Northing"), crs = 27700)
  
  pred_grid_inla<- read.csv("lsoa_prediction.csv")
  
  #overall average for LSOA prediction
  overall_avg <- pred_grid_inla %>%
    group_by(LSOA21CD, Easting, Northing) %>%
    summarise(pred_mean = mean(mean), .groups = "drop")
  
  overall_sd <- pred_grid_inla %>%
    group_by(LSOA21CD, Easting, Northing) %>%
    summarise(sd = mean(sd), .groups = "drop")
  
  

  lsoa_joined_pred_inla <- lsoa_boundaries %>%
    left_join(overall_avg, by = c("LSOA21CD" = "LSOA21CD")) %>%
    left_join(overall_sd, by = c("LSOA21CD" = "LSOA21CD")) %>%
    st_as_sf(coords = c("Easting", "Northing"), crs = 27700)
  
  # plot all together (INLA + RF)
  
  
  library(patchwork)
  library(ggplot2)
  
  # --- Plot 1: Predicted Mean ---
  p1 <- ggplot(lsoa_joined_pred_rf) +
    geom_sf(aes(fill = pred_mean), color = NA) +
    scale_fill_viridis_c(option = "plasma") +
    labs(title = "LSOA Prediction: RF (Mean)",
         fill = "NoV concentration\n(Log10gc per 100k capita)") +
    theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(),
          # --- Remove axis text and ticks ---
          axis.text.x       = element_blank(),
          axis.text.y       = element_blank(),
          axis.ticks        = element_blank())
  
  # --- Plot 2: Standard Deviation ---
  p2 <- ggplot(lsoa_joined_pred_rf) +
    geom_sf(aes(fill = sd), color = NA) +
    scale_fill_viridis_c(option = "mako") +
    labs(title = "LSOA Prediction: RF (SD)",
         fill = "Standard deviation") +
    theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(),
          # --- Remove axis text and ticks ---
          axis.text.x       = element_blank(),
          axis.text.y       = element_blank(),
          axis.ticks        = element_blank())
  
  # --- Combine & Save ---
  combined <- p1 + p2 +
    plot_layout(ncol = 2) +
    plot_annotation(title = "LSOA Random Forest Predictions")
  
  pdf("output/lsoa_prediction_plot_rf_combined.pdf", width = 16, height = 8)
  print(combined)
  dev.off()
  
}