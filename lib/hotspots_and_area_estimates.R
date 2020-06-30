
# Libraries ---------------------------------------------------------------

setwd("LandTrendr/drought/")

library(sf)
library(tidyverse)
library(raster)
library(patchwork)

# Prepare data ---------------------------------------------------------------

index <- "CWB"
thresholds <- read_csv(paste0("results/spatial_temporal/thresholds_", index, ".csv"))
threshold <- thresholds[thresholds$excess == "Excess mortality", ]$CWB_bp_mean
threshold_lower <- thresholds[thresholds$excess == "Excess mortality", ]$CWB_bp_lower
threshold_upper <- thresholds[thresholds$excess == "Excess mortality", ]$CWB_bp_upper

aic_temporal <- read_csv(paste0("results/spatial_temporal/aic_spatial_temporal_", index, ".csv"))

climate_residuals <- read_csv(paste0("data/", index, "_residuals.csv"))
climate_residuals <- climate_residuals %>% mutate_at(.vars = vars(lag, month), as.integer)

climate_residuals_select <- climate_residuals %>%
  filter(month == as.integer(aic_temporal[which.min(aic_temporal$aic), "month"]) & 
           lag == as.integer(aic_temporal[which.min(aic_temporal$aic), "lag"]))

climate_residuals_select <- climate_residuals_select %>%
  group_by(climate_id) %>%
  mutate(index_scaled = (index - mean(index)) / sd(index)) %>%
  ungroup()

dat <- read_csv("temp/dat.csv")

hotspots <- climate_residuals_select %>% 
  left_join(dat %>% group_by(climate_id) %>% 
              summarize(forest_ha = mean(forest_ha)) %>% 
              ungroup(), 
            by = c("climate_id")) %>%
  mutate(excess = ifelse(residuals_scaled > 0, 1, 0),
         drought = ifelse(index_scaled < threshold, 1, 0)) %>% 
  group_by(year, climate_id, ecoregions_name) %>% 
  summarize(total_total_mortality_ha = sum(rate * forest_ha),
            expected_mortality_ha = sum((excess * drought) * (fitted * forest_ha)),
            total_mortality_ha = sum((excess * drought) * (rate * forest_ha)),
            excess_mortality_proportion = sum(excess * residuals_scaled * drought)) %>%
  ungroup() %>%
  mutate(excess_mortality_ha = total_mortality_ha - expected_mortality_ha)

hotspots_lower <- climate_residuals_select %>% 
  left_join(dat %>% group_by(climate_id) %>% 
              summarize(forest_ha = mean(forest_ha)) %>% 
              ungroup(), 
            by = c("climate_id")) %>%
  mutate(excess = ifelse(residuals_scaled > 0, 1, 0),
         drought = ifelse(index_scaled < threshold_lower, 1, 0)) %>% 
  group_by(year, climate_id, ecoregions_name) %>% 
  summarize(total_total_mortality_ha = sum(rate * forest_ha),
            expected_mortality_ha = sum((excess * drought) * (fitted * forest_ha)),
            total_mortality_ha = sum((excess * drought) * (rate * forest_ha)),
            excess_mortality_proportion = sum(excess * residuals_scaled * drought)) %>%
  ungroup() %>%
  mutate(excess_mortality_ha = total_mortality_ha - expected_mortality_ha)

hotspots_upper <- climate_residuals_select %>% 
  left_join(dat %>% group_by(climate_id) %>% 
              summarize(forest_ha = mean(forest_ha)) %>% 
              ungroup(), 
            by = c("climate_id")) %>%
  mutate(excess = ifelse(residuals_scaled > 0, 1, 0),
         drought = ifelse(index_scaled < threshold_upper, 1, 0)) %>% 
  group_by(year, climate_id, ecoregions_name) %>% 
  summarize(total_total_mortality_ha = sum(rate * forest_ha),
            expected_mortality_ha = sum((excess * drought) * (fitted * forest_ha)),
            total_mortality_ha = sum((excess * drought) * (rate * forest_ha)),
            excess_mortality_proportion = sum(excess * residuals_scaled * drought)) %>%
  ungroup() %>%
  mutate(excess_mortality_ha = total_mortality_ha - expected_mortality_ha)

# Mapping hotspots --------------------------------------------------------

countries <- read_sf("data/countries_europe.shp")

reference_grid_shp <- read_sf("data/reference_grid.shp")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, projection(countries))
st_crs(world) <- st_crs(countries)
world <- st_crop(world, st_bbox(countries) + c(-0.05, -0.05, 0.01, 0.01) * as.double(st_bbox(countries)))

# Annual

hotspots_grid <- reference_grid_shp %>%
  right_join(hotspots, by = "climate_id") %>%
  mutate(excess_mortality_proportion_max = ifelse(excess_mortality_proportion >= 1, 1, excess_mortality_proportion))

p_hotspots <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = hotspots_grid, 
          mapping = aes(fill = excess_mortality_proportion_max * 100), color = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  theme_linedraw() +
  scale_fill_gradient(low = "white", high = "#b2182b", 
                      breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", ">100")) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_colorbar(title = "Drought related excess mortality (%)\n", 
                               title.position = "bottom",
                               title.align = 0.5))  +
  coord_sf(expand = FALSE) +
  facet_wrap(~year) +
  geom_text(data = data.frame(year = 1987:2016), aes(x = Inf, y = Inf, label = year),
            hjust = 1.25, vjust = 5, size = 2.5)

# Lower bound

hotspots_grid_lower <- reference_grid_shp %>%
  right_join(hotspots_lower, by = "climate_id") %>%
  mutate(excess_mortality_proportion_max = ifelse(excess_mortality_proportion >= 1, 1, excess_mortality_proportion))

p_hotspots_lower <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = hotspots_grid_lower, 
          mapping = aes(fill = excess_mortality_proportion_max * 100), color = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  theme_linedraw() +
  scale_fill_gradient(low = "white", high = "#b2182b", 
                      breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", ">100")) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_colorbar(title = "Drought related excess mortality (%)\n", 
                               title.position = "bottom",
                               title.align = 0.5))  +
  coord_sf(expand = FALSE) +
  facet_wrap(~year) +
  geom_text(data = data.frame(year = 1987:2016), aes(x = Inf, y = Inf, label = year),
            hjust = 1.25, vjust = 5, size = 2.5)

# Upper bound

hotspots_grid_upper <- reference_grid_shp %>%
  right_join(hotspots_upper, by = "climate_id") %>%
  mutate(excess_mortality_proportion_max = ifelse(excess_mortality_proportion >= 1, 1, excess_mortality_proportion))

p_hotspots_upper <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = hotspots_grid_upper, 
          mapping = aes(fill = excess_mortality_proportion_max * 100), color = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  theme_linedraw() +
  scale_fill_gradient(low = "white", high = "#b2182b", 
                      breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", ">100")) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_colorbar(title = "Drought related excess mortality (%)\n", 
                               title.position = "bottom",
                               title.align = 0.5))  +
  coord_sf(expand = FALSE) +
  facet_wrap(~year) +
  geom_text(data = data.frame(year = 1987:2016), aes(x = Inf, y = Inf, label = year),
            hjust = 1.25, vjust = 5, size = 2.5)

# Save all maps

ggsave("results/hotspots.pdf", p_hotspots, width = 7.5, height = 7.5)
ggsave("results/hotspots_lower.pdf", p_hotspots_lower, width = 7.5, height = 7.5)
ggsave("results/hotspots_upper.pdf", p_hotspots_upper, width = 7.5, height = 7.5)

# Cummulativ

hotspots_grid_cummulative <- reference_grid_shp %>%
  right_join(hotspots %>%
               group_by(climate_id) %>%
               summarize(excess_mortality_ha = sum(excess_mortality_ha, na.rm = TRUE),
                         total_total_mortality_ha = sum(total_total_mortality_ha, na.rm = TRUE)), 
             by = "climate_id") %>%
  mutate(excess_mortality_proportion = excess_mortality_ha / total_total_mortality_ha) %>%
  mutate(excess_mortality_proportion_max = ifelse(excess_mortality_proportion >= 0.3, 0.3, excess_mortality_proportion))

hotspots_grid_cummulative_upper <- reference_grid_shp %>%
  right_join(hotspots_upper %>%
               group_by(climate_id) %>%
               summarize(excess_mortality_ha = sum(excess_mortality_ha, na.rm = TRUE),
                         total_total_mortality_ha = sum(total_total_mortality_ha, na.rm = TRUE)), 
             by = "climate_id") %>%
  mutate(excess_mortality_proportion = excess_mortality_ha / total_total_mortality_ha) %>%
  mutate(excess_mortality_proportion_max = ifelse(excess_mortality_proportion >= 0.3, 0.3, excess_mortality_proportion))

hotspots_grid_cummulative_lower <- reference_grid_shp %>%
  right_join(hotspots_lower %>%
               group_by(climate_id) %>%
               summarize(excess_mortality_ha = sum(excess_mortality_ha, na.rm = TRUE),
                         total_total_mortality_ha = sum(total_total_mortality_ha, na.rm = TRUE)), 
             by = "climate_id") %>%
  mutate(excess_mortality_proportion = excess_mortality_ha / total_total_mortality_ha) %>%
  mutate(excess_mortality_proportion_max = ifelse(excess_mortality_proportion >= 0.3, 0.3, excess_mortality_proportion))

p <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = hotspots_grid_cummulative, 
          mapping = aes(fill = excess_mortality_proportion_max * 100), color = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  theme_linedraw() +
  scale_fill_gradient(low = "white", high = "#b2182b", labels = c(c(seq(0, 20, 10), "> 30"))) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
        plot.title = element_text(size = 9)) +
  labs(title = "Percent of total mortality attributable to drought") +
  coord_sf(expand = FALSE)

p_lower <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = hotspots_grid_cummulative_lower, 
          mapping = aes(fill = excess_mortality_proportion_max * 100), color = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  theme_linedraw() +
  scale_fill_gradient(low = "white", high = "#b2182b", labels = c(c(seq(0, 20, 10), "> 30"))) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
        plot.title = element_text(size = 9)) +
  labs(title = "Percent of total mortality attributable to drought") +
  coord_sf(expand = FALSE)

p_upper <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = hotspots_grid_cummulative_upper, 
          mapping = aes(fill = excess_mortality_proportion_max * 100), color = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  theme_linedraw() +
  scale_fill_gradient(low = "white", high = "#b2182b", labels = c(c(seq(0, 20, 10), "> 30"))) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
        plot.title = element_text(size = 9)) +
  labs(title = "Percent of total mortality attributable to drought") +
  coord_sf(expand = FALSE)

ggsave("results/hotspots_cummulative.pdf", p, width = 5.5, height = 5.5)
ggsave("results/hotspots_cummulative_lower.pdf", p_lower, width = 5.5, height = 5.5)
ggsave("results/hotspots_cummulative_upper.pdf", p_upper, width = 5.5, height = 5.5)

hotspots_grid %>%
  st_drop_geometry() %>%
  summarise(excess_mortality_percent = mean(excess_mortality_proportion) * 100)

hotspots %>%
  summarize(excess_mortality_ha = sum(excess_mortality_ha, na.rm = TRUE),
            total_total_mortality_ha = sum(total_total_mortality_ha, na.rm = TRUE)) %>%
  mutate(excess_mortality_proportion = excess_mortality_ha / total_total_mortality_ha * 100)

hotspots_lower %>%
  summarize(excess_mortality_ha = sum(excess_mortality_ha, na.rm = TRUE),
            total_total_mortality_ha = sum(total_total_mortality_ha, na.rm = TRUE)) %>%
  mutate(excess_mortality_proportion = excess_mortality_ha / total_total_mortality_ha * 100)

hotspots_upper %>%
  summarize(excess_mortality_ha = sum(excess_mortality_ha, na.rm = TRUE),
            total_total_mortality_ha = sum(total_total_mortality_ha, na.rm = TRUE)) %>%
  mutate(excess_mortality_proportion = excess_mortality_ha / total_total_mortality_ha * 100)
  
# Area estimates ----------------------------------------------------------

hotspots_tmp <- vector("list", length = 3)

k <- 0

thresnames <- c("Mean", "Lower", "Upper")

for (threshold in thresholds[1, c(2, 3, 4)]) {
  
  k <- k + 1
  
  hotspots_tmp[[k]] <- climate_residuals_select %>% 
    left_join(dat %>% group_by(climate_id) %>% 
                summarize(forest_ha = mean(forest_ha)) %>% 
                ungroup(), 
              by = c("climate_id")) %>%
    mutate(excess = ifelse(residuals_scaled > 0, 1, 0),
           drought = ifelse(index_scaled < threshold, 1, 0)) %>% 
    group_by(year, climate_id, ecoregions_name) %>% 
    summarize(expected_mortality_ha = sum((excess * drought) * (fitted * forest_ha)),
              total_mortality_ha = sum((excess * drought) * (rate * forest_ha)),
              excess_mortality_proportion = sum(excess * residuals_scaled * drought)) %>%
    ungroup() %>%
    mutate(excess_mortality_ha = total_mortality_ha - expected_mortality_ha) %>%
    mutate(threshold = thresnames[k])
  
}

hotspots_uncertainty <- hotspots_tmp %>%
  bind_rows()

hotspots_uncertainty <- hotspots_uncertainty %>%
  group_by(year, threshold) %>%
  summarize(excess_mortality_ha = sum(excess_mortality_ha)) %>%
  spread(key = threshold, value = excess_mortality_ha)

p_areaaffected <- hotspots_uncertainty %>%
  ggplot(.) + 
  geom_bar(aes(x = year, y = Mean), stat = "identity", fill = "#b2182b") + 
  geom_errorbar(aes(x = year, ymin = Lower, ymax = Upper), width = 0.25) + 
  theme_linedraw() +
  theme(legend.position = "none",
        legend.justification = c(0, 1),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", hjust = 0, size = 9),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.4), "cm")) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(labels = scales::comma, expand = c(0, 1), limits = c(0, 180000)) +
  labs(x = "Year", y = "Hectares per year", col = NULL, title = "Drought-related excess mortality")

ggsave("results/area_affected.pdf", p_areaaffected, width = 3.5, height = 3.5)

hotspots_uncertainty %>%
  gather(key = uncertainty, value = excess_mortality_ha, -year) %>%
  group_by(uncertainty) %>%
  summarize(sum_ha = sum(excess_mortality_ha)) %>%
  mutate(percent = sum_ha / sum(dat$disturbance_ha) * 100)

hotspots_uncertainty %>%
  gather(key = uncertainty, value = excess_mortality_ha, -year) %>%
  group_by(year, uncertainty) %>%
  summarize(excess_mortality_ha = sum(excess_mortality_ha)) %>%
  spread(key = uncertainty, value = excess_mortality_ha) %>%
  View(.)

# Spatial analysis --------------------------------------------------------

reference_grid_shp <- read_sf("data/reference_grid.shp")
countries <- read_sf("data/countries_europe.shp")
st_crs(countries) <- st_crs(reference_grid_shp)

reference_grid_shp <- st_crop(reference_grid_shp, countries)
reference_grid_shp_countries <- st_intersection(reference_grid_shp, countries)

strata <- read_delim("../TimeSync/data/strata/country_grouping.csv", delim = ";")

strata <- reference_grid_shp_countries %>%
  st_drop_geometry() %>%
  left_join(strata, by = c("ISO_CC" = "iso_code")) %>%
  dplyr::select(climate_id, country_name, euro_region)

hotspots_tmp <- vector("list", length = 3)

k <- 0

thresnames <- c("Mean", "Lower", "Upper")

for (threshold in thresholds[1, c(2, 3, 4)]) {
  
  k <- k + 1
  
  hotspots_tmp[[k]] <- climate_residuals_select %>% 
    left_join(dat %>% group_by(climate_id) %>% 
                summarize(forest_ha = mean(forest_ha)) %>% 
                ungroup(), 
              by = c("climate_id")) %>%
    left_join(strata, by = "climate_id") %>%
    mutate(excess = ifelse(residuals_scaled > 0, 1, 0),
           drought = ifelse(index_scaled < threshold, 1, 0)) %>% 
    group_by(year, climate_id, euro_region) %>% 
    summarize(expected_mortality_ha = sum((excess * drought) * (fitted * forest_ha)),
              total_mortality_ha = sum((excess * drought) * (rate * forest_ha)),
              excess_mortality_proportion = sum(excess * residuals_scaled * drought)) %>%
    ungroup() %>%
    mutate(excess_mortality_ha = total_mortality_ha - expected_mortality_ha) %>%
    mutate(threshold = thresnames[k])
  
}

hotspots_tmp %>%
  bind_rows() %>%
  group_by(threshold, euro_region) %>%
  summarize(excess_mortality_ha = sum(excess_mortality_ha)) %>%
  spread(key = threshold, value = excess_mortality_ha) %>%
  ungroup() %>%
  mutate(Lower = Lower / sum(Lower),
         Mean = Mean / sum(Mean),
         Upper = Upper / sum(Upper))

# Sensitivity analysis ----------------------------------------------------

hotspots_tmp <- vector("list")

k <- 0

for (threshold in c(-0.5, -1, -1.5, -2)) {
  
  k <- k + 1
  
  hotspots_tmp[[k]] <- climate_residuals_select %>% 
    left_join(dat %>% group_by(climate_id) %>% 
                summarize(forest_ha = mean(forest_ha)) %>% 
                ungroup(), 
              by = c("climate_id")) %>%
    mutate(excess = ifelse(residuals_scaled > 0, 1, 0),
           drought = ifelse(index_scaled < threshold, 1, 0)) %>% 
    group_by(year, climate_id, ecoregions_name) %>% 
    summarize(expected_mortality_ha = sum((excess * drought) * (fitted * forest_ha)),
              total_mortality_ha = sum((excess * drought) * (rate * forest_ha)),
              excess_mortality_proportion = sum(excess * residuals_scaled * drought)) %>%
    ungroup() %>%
    mutate(excess_mortality_ha = total_mortality_ha - expected_mortality_ha) %>%
    mutate(threshold = paste0("CWB < ", threshold, " SD"))
  
}

hotspots_sensitivity <- hotspots_tmp %>%
  bind_rows()

p_sensitivity <- hotspots_sensitivity %>%
  group_by(year, threshold) %>%
  summarize(excess_mortality_ha = sum(excess_mortality_ha)) %>%
  ggplot(., aes(x = year, y = excess_mortality_ha)) + 
  geom_bar(stat = "identity", fill = "#b2182b") + 
  theme_linedraw() +
  theme(legend.position = "none",
        legend.justification = c(1, 1),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", hjust = 0, size = 9),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.4), "cm")) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(labels = scales::comma, expand = c(0, 1)) +
  labs(x = "Year", y = "Hectares per year", col = NULL, title = "Drought-related excess mortality") +
  facet_wrap(~threshold, ncol = 4)

ggsave("results/area_affected_sensitivity.pdf", p_sensitivity, width = 8, height = 2.75)


# Test thresholds

threshold <- -1.27

climate_residuals_select %>% 
  left_join(dat %>% group_by(climate_id) %>% 
              summarize(forest_ha = mean(forest_ha)) %>% 
              ungroup(), 
            by = c("climate_id")) %>%
  mutate(excess = ifelse(residuals_scaled > 0, 1, 0),
         drought = ifelse(index_scaled < threshold, 1, 0)) %>% 
  group_by(year, climate_id, ecoregions_name) %>% 
  summarize(expected_mortality_ha = sum((excess * drought) * (fitted * forest_ha)),
            total_mortality_ha = sum((excess * drought) * (rate * forest_ha)),
            excess_mortality_proportion = sum(excess * residuals_scaled * drought)) %>%
  ungroup() %>%
  mutate(excess_mortality_ha = total_mortality_ha - expected_mortality_ha) %>%
  mutate(threshold = paste0("CWB < ", threshold, " SD")) %>%
  group_by(year, threshold) %>%
  summarize(excess_mortality_ha = sum(excess_mortality_ha)) %>%
  ggplot(., aes(x = year, y = excess_mortality_ha)) + 
  geom_bar(stat = "identity", fill = "#b2182b") + 
  theme_linedraw() +
  theme(legend.position = "none",
        legend.justification = c(1, 1),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", hjust = 0, size = 9),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.4), "cm")) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(labels = scales::comma, expand = c(0, 1)) +
  labs(x = "Year", y = "Hectares per year", col = NULL, title = "Drought-related excess mortality")

