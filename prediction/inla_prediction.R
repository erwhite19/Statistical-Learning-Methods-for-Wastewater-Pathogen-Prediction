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



inla.setOption(num.threads= "6")

## read data -------------------------------------------------------------------
nov_df <- read.csv("nov_df_reduced.csv")
pred_grid <- fread("data/lsoa_grid_prediction.csv")

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
posterior.samples = inla.posterior.sample(n = 1000, result=pred.lsoa, intern = FALSE, use.improved.mean = TRUE, add.names = TRUE, seed = 0L)

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
    lower   = pred.lsoa$summary.fitted.values[index, "0.025quant"],
    upper   = pred.lsoa$summary.fitted.values[index, "0.975quant"],
    sd        = pred.lsoa$summary.fitted.values[index, "sd"]
  )


# Create the column names for columns 1 to 100
new_names <- paste0("sample_", 1:1000)

# Assign the new names to columns 1 to 100
colnames(samples.fitted)[1:1000] <- new_names
columns_to_convert <- 1:1000

# Step 3: Convert specified columns to numeric
pred_grid<- as.data.frame(samples.fitted)
pred_grid[, columns_to_convert] <- apply(pred_grid[, columns_to_convert], 2, as.numeric)

fwrite(samples.df, "output/inla_predictions.csv")

pred_grid<- fread("output/inla_predictions.csv")

# lsoa population
lsoa_pop<- read.csv("Data/LSOA/population_lsoa_21.csv")

# Convert both to data.table (does this in place, no copy made)
setDT(pred_grid)
setDT(lsoa_pop)

sample_cols <- grep("^sample_", names(pred_grid), value = TRUE)


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


write.csv(weekly_summary, "output/inla_nat_pred.csv")

# Get LSOA mean (use data.table to help with large dataset)


lsoa_samples <- pred_grid[  # left join
  ,                                                         # no row filter
  lapply(.SD, function(x) mean(x, w = population, na.rm = TRUE)),
  by = LSOA21CD,                                                # group by week
  .SDcols = sample_cols                                     # only sample columns
]


# compute intervals
sample_cols <- grep("^sample_", names(lsoa_samples), value = TRUE)
sample_mat <- as.matrix(lsoa_samples[, ..sample_cols])

lsoa_summary <- data.frame(
  LSOA21CD     = lsoa_samples$LSOA21CD,
  mean = rowMeans(sample_mat),
  lower  = apply(sample_mat, 1, quantile, 0.025),
  upper  = apply(sample_mat, 1, quantile, 0.975)
)


write.csv(lsoa_summary, "output/inla_lsoa_pred.csv")

