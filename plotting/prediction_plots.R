# =============================================================================
# Norovirus Prediction Visualisation - Spatial Maps and National Trends
# =============================================================================

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(sf)

# -----------------------------------------------------------------------------
# Load Prediction Data
# -----------------------------------------------------------------------------

# Random Forest predictions
pred_rf      <- fread("outputs/rf_prediction_mean.csv")
nat_pred_rf  <- read.csv("outputs/rf_nat_pred.csv")
space_pred_rf <- read.csv("outputs/rf_lsoa_prediction.csv")

# INLA predictions
pred_inla      <- fread("outputs/inla_prediction_mean.csv")
nat_pred_inla  <- read.csv("outputs/inla_nat_pred.csv")
space_pred_inla <- read.csv("outputs/inla_lsoa_prediction.csv")

# LSOA boundaries (England only)
lsoa_boundaries <- st_read("data/LSOA2021_boundaries/LSOA2021_boundaries.shp") %>%
  dplyr::select(LSOA21CD, geometry) %>%
  dplyr::filter(grepl("^E0", LSOA21CD))

# LSOA population
lsoa_pop <- read.csv("data/population_lsoa_21.csv")

# -----------------------------------------------------------------------------
# Join Predictions to Boundaries
# -----------------------------------------------------------------------------

space_rf_sf <- lsoa_boundaries %>%
  left_join(space_pred_rf, by = "LSOA21CD")

space_inla_sf <- lsoa_boundaries %>%
  left_join(space_pred_inla, by = "LSOA21CD")

# -----------------------------------------------------------------------------
# Shared Plot Theme
# -----------------------------------------------------------------------------

map_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks       = element_blank()
  )

# -----------------------------------------------------------------------------
# LSOA Spatial Maps - Mean and SD
# -----------------------------------------------------------------------------

# Shared fill limits for consistent legends across models
global_min_mean <- min(space_inla_sf$mean, space_rf_sf$mean, na.rm = TRUE)
global_max_mean <- max(space_inla_sf$mean, space_rf_sf$mean, na.rm = TRUE)
global_min_sd   <- min(space_inla_sf$sd,   space_rf_sf$sd,   na.rm = TRUE)
global_max_sd   <- max(space_inla_sf$sd,   space_rf_sf$sd,   na.rm = TRUE)

nov_conc_label <- expression(atop(
  "Norovirus concentration",
  (log[10] ~ gene ~ copies ~ "/ 100k capita")
))

# INLA: predicted mean
p1 <- ggplot(space_inla_sf) +
  geom_sf(aes(fill = mean), color = NA, linewidth = 0) +
  scale_fill_viridis_c(option = "plasma",
                       limits = c(global_min_mean, global_max_mean)) +
  labs(fill = nov_conc_label) +
  map_theme

# INLA: standard deviation
p2 <- ggplot(space_inla_sf) +
  geom_sf(aes(fill = sd), color = NA, linewidth = 0) +
  scale_fill_viridis_c(option = "mako",
                       limits = c(global_min_sd, global_max_sd)) +
  labs(fill = "Standard deviation") +
  map_theme

# RF: predicted mean
p3 <- ggplot(space_rf_sf) +
  geom_sf(aes(fill = mean), color = NA, linewidth = 0) +
  scale_fill_viridis_c(option = "plasma",
                       limits = c(global_min_mean, global_max_mean)) +
  labs(fill = nov_conc_label) +
  map_theme

# RF: standard deviation
p4 <- ggplot(space_rf_sf) +
  geom_sf(aes(fill = sd), color = NA, linewidth = 0) +
  scale_fill_viridis_c(option = "mako",
                       limits = c(global_min_sd, global_max_sd)) +
  labs(fill = "Standard deviation") +
  map_theme

# Legend formatting for bottom row
legend_guide <- guides(fill = guide_colourbar(
  barwidth       = 10,
  barheight      = 0.8,
  title.position = "top",
  title.hjust    = 0.5
))

legend_theme <- theme(
  legend.position  = "bottom",
  legend.direction = "horizontal"
)

combined <- (
  (p1 + theme(legend.position = "none")) +
    (p2 + theme(legend.position = "none"))
) / (
  (p3 + legend_guide + legend_theme) +
    (p4 + legend_guide + legend_theme)
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

png("outputs/plots/lsoa_maps.png",
    width = 10, height = 12, units = "in", res = 300, type = "cairo")
print(combined)
dev.off()

# -----------------------------------------------------------------------------
# National Trend Plot
# -----------------------------------------------------------------------------

nat_pred_inla <- nat_pred_inla %>%
  mutate(date = start_date + 7 * (week - 1))

nat_pred_rf <- nat_pred_rf %>%
  mutate(date = start_date + 7 * (week - 1))

nat_trends <- ggplot() +
  
  # Credible / confidence ribbons
  geom_ribbon(data = nat_pred_inla,
              aes(x = date, ymin = lower, ymax = upper),
              fill = "#f4a582", alpha = 0.4) +
  geom_ribbon(data = nat_pred_rf,
              aes(x = date, ymin = lower, ymax = upper),
              fill = "#92a8d1", alpha = 0.4) +
  
  # Mean lines
  geom_line(data = nat_pred_inla,
            aes(x = date, y = mean, color = "INLA"),
            linewidth = 0.8) +
  geom_line(data = nat_pred_rf,
            aes(x = date, y = mean, color = "RF"),
            linewidth = 1.0) +
  
  # Colour scale
  scale_color_manual(
    name   = NULL,
    values = c("INLA" = "#d6604d", "RF" = "#4575b4"),
    guide  = guide_legend(override.aes = list(linewidth = 1.2))
  ) +
  
  # Axes
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand      = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = seq(3.8, 6.0, by = 0.2),
    limits = c(3.6, 6.1)
  ) +
  
  # Labels
  labs(
    x = "Date",
    y = expression(atop(
      "Norovirus concentration",
      log[10] ~ "(gene copies / 100k capita)"
    ))
  ) +
  
  # Theme
  theme_classic(base_size = 13) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.position   = c(0.08, 0.92),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.text       = element_text(size = 12),
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black")
  )

png("outputs/plots/nat_trends.png",
    width = 12, height = 6, units = "in", res = 300)
print(nat_trends)
dev.off()

# -----------------------------------------------------------------------------
# Spatio-Temporal Maps - Selected Weeks
# -----------------------------------------------------------------------------

start_date    <- as.Date("2021-05-24")
weeks_to_plot <- seq(1, 46, by = 8)

# Add date column to spatial data
lsoa_inla <- lsoa_boundaries %>%
  left_join(pred_inla, by = "LSOA21CD") %>%
  mutate(date = start_date + 7 * (week - 1))

lsoa_rf <- lsoa_boundaries %>%
  left_join(pred_rf, by = "LSOA21CD") %>%
  mutate(date = start_date + 7 * (week - 1))

fill_limits <- c(
  floor(min(c(lsoa_inla$mean,  lsoa_rf$mean),  na.rm = TRUE)),
  ceiling(max(c(lsoa_inla$mean, lsoa_rf$mean), na.rm = TRUE))
)

st_map_theme <- theme_minimal() +
  theme(
    legend.position  = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    plot.title       = element_text(size = 8, hjust = 0.5),
    plot.margin      = margin(1, 1, 1, 1)
  )

plots_inla <- list()
plots_rf   <- list()

for (t in weeks_to_plot) {
  
  week_data_inla <- subset(lsoa_inla, week == t)
  week_data_rf   <- subset(lsoa_rf,   week == t)
  date_label     <- format(unique(week_data_inla$date)[1], "%d %b %Y")
  
  plots_inla[[paste0("week_", t)]] <- ggplot() +
    geom_sf(data = week_data_inla, aes(fill = mean), color = NA) +
    scale_fill_viridis_c(option = "plasma", limits = fill_limits) +
    labs(title = date_label) +
    st_map_theme
  
  plots_rf[[paste0("week_", t)]] <- ggplot() +
    geom_sf(data = week_data_rf, aes(fill = mean), color = NA) +
    scale_fill_viridis_c(option = "plasma", limits = fill_limits) +
    labs(title = date_label) +
    st_map_theme
}

# Combine into labelled panels
panel_inla <- wrap_plots(plots_inla, ncol = 3, nrow = 2) +
  plot_annotation(
    title = "A",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0))
  )

panel_rf <- wrap_plots(plots_rf, ncol = 3, nrow = 2) +
  plot_annotation(
    title = "B",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0))
  )

# Shared legend
legend_dummy <- ggplot(
  data.frame(x = seq(fill_limits[1], fill_limits[2], length.out = 100)),
  aes(x = x, y = 1, fill = x)
) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", limits = fill_limits) +
  labs(fill = nov_conc_label) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 9, hjust = 0.5),
    legend.text     = element_text(size = 8)
  )

shared_legend <- cowplot::get_legend(legend_dummy)

st_plot <- plot_grid(
  panel_inla,
  panel_rf,
  shared_legend,
  ncol        = 1,
  rel_heights = c(3, 3, 0.3)
)

png("outputs/plots/st_maps.png",
    width = 10, height = 12, units = "in", res = 300, type = "cairo")
print(st_plot)
dev.off()

