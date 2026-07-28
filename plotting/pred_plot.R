


# RF data

pred_rf<= fread<- (".../rf_lsoa_prediction.csv")

# aggregate for plotting

# national trend

# spatial mean

# space time


pred_inla<- read.csv(" .../inla_lsoa_prediction.csv")

# aggregate for plotting

# national trend

# spatial mean

# space time

  
  
# National trends





# spatial mean
pred_grid_rf<- fread("output/rf_lsoa_prediction_new.csv")

# add population



#overall average for LSOA prediction
overall_avg_rf <- pred_grid_rf %>%
  group_by(LSOA21CD, Easting, Northing) %>%
  summarise(pred_mean = mean(predicted), .groups = "drop")

overall_sd_rf <- pred_grid_rf %>%
  group_by(LSOA21CD, Easting, Northing) %>%
  summarise(sd = mean(sd), .groups = "drop")



lsoa_joined_pred_rf <- lsoa_boundaries %>%
  left_join(overall_avg_rf, by = c("LSOA21CD" = "LSOA21CD")) %>%
  left_join(overall_sd_rf, by = c("LSOA21CD" = "LSOA21CD")) %>%
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700)

pred_grid_inla<- fread("lsoa_prediction.csv")

#overall average for LSOA prediction
overall_avg_inla <- pred_grid_inla %>%
  group_by(LSOA21CD, Easting, Northing) %>%
  summarise(pred_mean = mean(pred_mean), .groups = "drop")

overall_sd_inla <- pred_grid_inla %>%
  group_by(LSOA21CD, Easting, Northing) %>%
  summarise(sd = mean(sd), .groups = "drop")

# global min and max
global_min_mean<- min(overall_avg_inla$pred_mean, overall_avg_rf$pred_mean)
global_max_mean<- max(overall_avg_inla$pred_mean, overall_avg_rf$pred_mean)

global_min_sd<- min(overall_sd_inla$sd, overall_sd_rf$sd)
global_max_sd<- max(overall_sd_inla$sd, overall_sd_rf$sd)


lsoa_joined_pred_inla <- lsoa_boundaries %>%
  left_join(overall_avg_inla, by = c("LSOA21CD" = "LSOA21CD")) %>%
  left_join(overall_sd_inla, by = c("LSOA21CD" = "LSOA21CD")) %>%
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700)

# --- Plot 1: Predicted Mean ---
p1 <- ggplot(lsoa_joined_pred_inla) +
  geom_sf(aes(fill = pred_mean), color = NA, lwd = 0, linewidth = 0) +
  scale_fill_viridis_c(option = "plasma", 
                       limits = c(global_min_mean, global_max_mean)) +
  labs(fill = expression(atop("NoV concentration",
                              (log[10]~gene~copies~"/ 100k capita")))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank())

# --- Plot 2: Standard Deviation ---
p2 <- ggplot(lsoa_joined_pred_inla) +
  geom_sf(aes(fill = sd), color = NA, lwd = 0, linewidth = 0) +
  scale_fill_viridis_c(option = "mako", limits = c(global_min_sd, global_max_sd)) +
  labs(#title = "LSOA Prediction: INLA  +
    fill = "Standard deviation") +
  theme_bw() +
  theme(panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        # --- Remove axis text and ticks ---
        axis.text.x       = element_blank(),
        axis.text.y       = element_blank(),
        axis.ticks        = element_blank())
# --- Top row: NO legends ---
p1_no_leg <- p1 + theme(legend.position = "none")
p1_no_leg <- p2 + theme(legend.position = "none")

# --- Plot 1: Predicted Mean ---
p3 <- ggplot(lsoa_joined_pred_rf) +
  geom_sf(aes(fill = pred_mean), color = NA, lwd = 0, linewidth = 0) +
  scale_fill_viridis_c(option = "plasma", limits = c(global_min_mean, global_max_mean)) +
  labs(fill = expression(atop("NoV concentration",
                              (log[10]~gene~copies~"/ 100k capita"))))  +
  theme_bw() +
  theme(panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        # --- Remove axis text and ticks ---
        axis.text.x       = element_blank(),
        axis.text.y       = element_blank(),
        axis.ticks        = element_blank())

# --- Plot 2: Standard Deviation ---
p4 <- ggplot(lsoa_joined_pred_rf) +
  geom_sf(aes(fill = sd), color = NA, lwd = 0, linewidth = 0) +
  scale_fill_viridis_c(option = "mako", limits = c(global_min_sd, global_max_sd)) +
  labs(#title = "LSOA Prediction: RF (SD)",
    fill = "Standard deviation") +
  theme_bw() +
  theme(panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        # --- Remove axis text and ticks ---
        axis.text.x       = element_blank(),
        axis.text.y       = element_blank(),
        axis.ticks        = element_blank())





# --- Bottom row: legends centered below each plot ---
p3_leg <- p3 + guides(fill = guide_colourbar(barwidth       = 10,
                                             barheight      = 0.8,
                                             title.position = "top",
                                             title.hjust    = 0.5)) +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal")

p4_leg <- p4 + guides(fill = guide_colourbar(barwidth       = 10,
                                             barheight      = 0.8,
                                             title.position = "top",
                                             title.hjust    = 0.5)) +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal")

# --- Combine ---
combined <- (p1_no_leg + p2_no_leg) /
  (p3_leg    + p3_leg) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))


png(".../combined.png", width  = 10, height = 12, units= "in", res= 300, type = "cairo")
print(combined)
dev.off()





# space time






  





