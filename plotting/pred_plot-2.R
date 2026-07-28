library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(sf)

setwd("~/GitHub/Statistical-Learning-Methods-for-Wastewater-Pathogen-Prediction")

## load data (all data sets produced by prediction script)

# RF data 

pred_rf<- fread("outputs/rf_prediction_mean.csv")

# aggregate for plotting

# national trend

nat_pred_rf<- read.csv("outputs/rf_nat_pred.csv")

# spatial mean
space_pred_rf<- read.csv("outputs/rf_lsoa_prediction.csv")

pred_inla<- fread("outputs/inla_prediction_mean.csv") 

#pred_inla<- fread("outputs/inla_lsoa_prediction_samples.csv") # run prediction script


# national trend

nat_pred_inla<- read.csv("outputs/inla_nat_pred.csv")

# spatial mean
space_pred_inla<- read.csv("outputs/inla_lsoa_prediction.csv")



# LSOA boundaries
lsoa_boundaries <- st_read("data/LSOA2021_boundaries/LSOA2021_boundaries.shp") %>%
  dplyr::select(LSOA21CD, geometry)

lsoa_boundaries<- lsoa_boundaries[grepl("^E0", lsoa_boundaries$LSOA21CD), ]

# add population
lsoa_pop<- read.csv("data/population_lsoa_21.csv")


  
  
pred_rf_sf <- lsoa_boundaries %>%
  left_join(space_pred_rf, by = c("LSOA21CD" = "LSOA21CD")) %>%
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700)




#pred_grid_inla<- fread("outputs/inla_lsoa_prediction_samples.csv")

# add population

space_rf_sf <- lsoa_boundaries %>%
  left_join(space_pred_rf, by = c("LSOA21CD" = "LSOA21CD")) %>%
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700)

space_inla_sf <- lsoa_boundaries %>%
  left_join(space_pred_inla, by = c("LSOA21CD" = "LSOA21CD")) %>%
  st_as_sf(coords = c("Easting", "Northing"), crs = 27700)


# global min and max (for consistent legend)
global_min_mean<- min(space_inla_sf$pred_mean, space_rf_sf$pred_mean)
global_max_mean<- max(space_inla_sf$pred_mean, space_rf_sf$pred_mean)

global_min_sd<- min(space_inla_sf$sd, space_rf_sf$sd)
global_max_sd<- max(space_inla_sf$sd, space_rf_sf$sd)




# --- Plot 1: INLA Predicted Mean ---
p1 <- ggplot(space_inla_sf) +
  geom_sf(aes(fill = mean), color = NA, lwd = 0, linewidth = 0) +
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
p2 <- ggplot(space_inla_sf) +
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
p2_no_leg <- p2 + theme(legend.position = "none")

# --- Plot 1: Predicted Mean ---
p3 <- ggplot(space_rf_sf) +
  geom_sf(aes(fill = mean), color = NA, lwd = 0, linewidth = 0) +
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
p4 <- ggplot(space_rf_sf) +
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


png("outputs/space_map.png", width  = 10, height = 12, units= "in", res= 300, type = "cairo")
print(combined)
dev.off()



# space time

# ══════════════════════════════════════════════════════════════════════════════
# 1. Prep data — deduplicate and join to boundaries
# ══════════════════════════════════════════════════════════════════════════════

weeks_to_plot <- seq(1, 46, by = 8)   # Plot every 8 weeks

# ── INLA ──────────────────────────────────────────────────────────────────────


lsoa_inla <- lsoa_boundaries %>%
  left_join(pred_inla, by = "LSOA21CD")

# ── RF ────────────────────────────────────────────────────────────────────────


lsoa_rf <- lsoa_boundaries %>%
  left_join(pred_rf, by = "LSOA21CD")
# ══════════════════════════════════════════════════════════════════════════════
# 2. Shared fill limits
# ══════════════════════════════════════════════════════════════════════════════

fill_limits <- c(
  floor(min(c(lsoa_inla$mean,  lsoa_rf$mean),  na.rm = TRUE)),
  ceiling(max(c(lsoa_inla$mean, lsoa_rf$mean), na.rm = TRUE))
)

# ══════════════════════════════════════════════════════════════════════════════
# 3. Build plot lists
# ══════════════════════════════════════════════════════════════════════════════

plots_inla <- list()
plots_rf   <- list()

lsoa_inla$date <-  as.Date("2021-05-24") + 7*(lsoa_inla$week - 1)
lsoa_rf$date<- as.Date("2021-05-24") + 7*(lsoa_rf$week - 1)
for (t in weeks_to_plot) {
  
  week_data_inla <- subset(lsoa_inla, week == t)
  week_data_rf   <- subset(lsoa_rf,   week == t)
  date_label     <- format(as.Date(unique(lsoa_inla$date)[1]), "%d %b %Y")
  
  # ── INLA ────────────────────────────────────────────────────────────────────
  plots_inla[[paste0("week_", t)]] <- ggplot() +
    geom_sf(data = week_data_inla, aes(fill = mean), color = NA) +
    scale_fill_viridis_c(option = "plasma", 
                         limits = fill_limits) +
    labs(title = date_label,
         fill = expression(atop("NoV concentration",
                                (log[10]~gene~copies~"/ 100k capita")))) +
    theme_minimal() +
    theme(
      legend.position  = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      plot.title       = element_text(size = 8, hjust = 0.5),
      plot.margin      = margin(1, 1, 1, 1)
    )
  
  # ── RF ──────────────────────────────────────────────────────────────────────
  plots_rf[[paste0("week_", t)]] <- ggplot() +
    geom_sf(data = week_data_rf, aes(fill = mean), color = NA) +
    scale_fill_viridis_c(option = "plasma", 
                         limits = fill_limits) +
    labs(title = date_label,
         fill = expression(atop("NoV concentration",
                                (log[10]~gene~copies~"/ 100k capita"))))+ 
    theme_minimal() +
    theme(
      legend.position  = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      plot.title       = element_text(size = 8, hjust = 0.5),
      plot.margin      = margin(1, 1, 1, 1)
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Combine into 4 x 3 panels with A / B labels top-left
# ══════════════════════════════════════════════════════════════════════════════

panel_inla <- wrap_plots(plots_inla, ncol = 3, nrow = 2) +
  plot_annotation(
    title = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0)
    )
  )

panel_rf <- wrap_plots(plots_rf, ncol = 3, nrow = 2) +
  plot_annotation(
    title = "B",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0)
    )
  )

# ══════════════════════════════════════════════════════════════════════════════
# 5. Shared legend
# ══════════════════════════════════════════════════════════════════════════════

legend_dummy <- ggplot(
  data.frame(x = seq(fill_limits[1], fill_limits[2], length.out = 100)),
  aes(x = x, y = 1, fill = x)
) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", 
                       limits = c(global_min_mean, global_max_mean)) +
  labs(fill = expression(atop("NoV concentration",
                              (log[10]~gene~copies~"/ 100k capita")))) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 9,  hjust = 0.5),
    legend.text     = element_text(size = 8)
  )

shared_legend <- cowplot::get_legend(legend_dummy)

# ══════════════════════════════════════════════════════════════════════════════
# 6. Stack: A panel / B panel / legend
# ══════════════════════════════════════════════════════════════════════════════

st_plot <- plot_grid(
  panel_inla,
  panel_rf,
  shared_legend,
  ncol        = 1,
  rel_heights = c(3, 3, 0.3)   # 3 rows each panel, small legend strip
)

# ══════════════════════════════════════════════════════════════════════════════
# 7. Save
# ══════════════════════════════════════════════════════════════════════════════

png("st_combined.png", width  = 10, height = 12, units= "in", res= 300, type = "cairo")
print(st_plot)
dev.off()






  





