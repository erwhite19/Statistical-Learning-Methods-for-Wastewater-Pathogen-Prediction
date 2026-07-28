library(sf)
library(dplyr)
library(ggplot2)

setwd("~/Term3-project")

# RF grid prediction plot
grid_RF <- read.csv("outputs/prediction/rf_lsoa_prediction.csv")

#create date and data_label for plotting
grid_RF$date <-  as.Date("2021-05-24") + 7*(grid_RF$date_index - 1)
grid_RF$date <- format(grid_RF$date, "%d/%m/%Y")
grid_RF$date_label <- factor(grid_RF$date, levels = unique(grid_RF$date))



# shapefile
england<- st_read("Data/shapefiles/england/england_crop.shp")
england<- st_transform(england, crs= 27700)

# create limits for the outcome so the plots are in the same scale
global_min <- min(grid$pred_mean, grid_RF$predicted)
global_max <- max(grid$pred_mean, grid_RF$predicted)


i <- seq(1,5)
subset_grid_RF <- grid_RF %>%
  filter(date_index == i)

p <- ggplot(subset_grid_RF, aes(x = Easting, y = Northing, color = predicted)) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", limits = c(global_min, global_max)) +
  geom_sf(data = england, fill = NA, color ="grey", inherit.aes = FALSE)+
  coord_sf() +
  theme_bw() +
  facet_wrap(~date_label, ncol = 5)+
  labs(title = paste("Spatial Distribution of RF Predictions by Week"),
       color = "predicted Value") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p)

week1 <- grid_RF %>% 
  filter(date_index == 1)

p1 <- ggplot(week1, aes(x = Easting, y = Northing, color = predicted)) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", limits = c(global_min, global_max)) +
  geom_sf(data = england, fill = NA, color ="grey", inherit.aes = FALSE)+
  coord_sf() +
  theme_bw() +
  labs(title = paste("Spatial Distribution of RF Predictions Week 1"),
       color = "Predicted Value") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p1)
