
library(tidyverse)


# Definition mortality pulse ----------------------------------------------

dat <- read_csv("temp/dat.csv")

climate_id_exp <- sample(unique(dat$climate_id), 1)
# Used for paper: 31072
climate_id_exp <- 31072

exp_glm <- dat %>%
  filter(climate_id == climate_id_exp) %>%
  mutate(., 
         disturbance_n = disturbance_ha / 0.09,
         forest_n = forest_ha / 0.09) %>%
  glm(cbind(disturbance_n, forest_n) ~ year, data = ., family = binomial(link = "logit"))

dat_exp <- dat %>%
  filter(climate_id == climate_id_exp) %>%
  mutate(fitted = fitted(exp_glm),
         residuals = rate - fitted,
         residuals_relative = residuals / fitted) %>%
  dplyr::select(climate_id, year, rate, residuals, residuals_relative, fitted)

p_exp1 <- ggplot(dat_exp, aes(x = year, y = rate * 100)) +
  geom_line() +
  geom_line(aes(y = fitted * 100, col = "Trend"), size = 1.2, linetype = "dashed") +
  theme_linedraw() +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        plot.margin = unit(c(0.5, 0.1, 0, 0.5), "cm")) +
  labs(y = bquote("Canopy mortality (%"~yr^-1*")"), col = NULL) +
  scale_color_manual(values = "lightgrey") +
  scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010, 2015), expand = c(0, 0))

p_exp2 <- ggplot(dat_exp, aes(x = year, y = residuals_relative * 100, fill = residuals_relative)) +
  geom_bar(stat = "identity") +
  theme_linedraw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.margin = unit(c(0, 0.1, 0.5, 0.5), "cm")) +
  labs(y = "Deviation from trend (%)", x = "Year") +
  scale_fill_gradient2(low = "#4885C1", mid = "white", high = "#AE3A4E") +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010, 2015), expand = c(0, 0))

p_exp <- p_exp1 + p_exp2 + plot_layout(ncol = 1, height = c(0.5, 0.5)) 

p_exp

ggsave(paste0("results/residuals_example_", climate_id_exp, ".pdf"), p_exp, 
       width = 5, height = 5)

# Lag-month ----------------------------------------------------------------

dat_aic <- list.files("results", pattern = glob2rx("*aic_temporal_*csv"), full.names = TRUE) %>%
  map(read_csv) %>%
  set_names(c("CWB", "PREC", "VPD")) %>%
  bind_rows(.id = "droughtindex")

dat_aic <- dat_aic %>%
  group_by(droughtindex) %>%
  mutate(d_aic = aic - min(aic)) %>%
  mutate(rel_lik = exp(-0.5 * d_aic)) %>%
  mutate(weight = rel_lik / sum(rel_lik))
  
p_aic <- ggplot(dat_aic,
                aes(x = str_pad(lag, 2, side = "left", pad = "0"),
                    y = str_pad(month, 2, side = "left", pad = "0"),
                    fill = aic)) +
  geom_tile() +
  scale_fill_viridis_c(breaks = c(min(dat_aic$aic), max(dat_aic$aic)),
                       labels = c("Lowest", "Highest")) +
  theme_linedraw() +
  theme(panel.spacing = unit(0.1, "cm"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 7),
        legend.position = "right",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.margin=margin(0, 0, 0, 0),
        legend.box.margin=margin(-10, 2, -10, -5),
        plot.title = element_text(size = 8),
        plot.margin = unit(c(0.2, 0, 0.2, 0.2), "cm")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Temporal averaging (months)", y = "Observation month", fill = "AIC") +
  guides(fill =  guide_colorbar(barheight = unit(1.75, "cm"),
                                barwidth = unit(0.15, "cm"), 
                                title.position = "top")) +
  facet_wrap(~droughtindex, scales = "free")

ggsave("figures/fig01.pdf", p_aic, width = 3.5, height = 1.55)


# Pulses ------------------------------------------------------------------

climate_residuals <- read.csv(paste0("data/CWB_residuals.csv"), stringsAsFactors = FALSE)

residuals <- climate_residuals %>%
  filter(lag == 1 & month == 1) %>%
  dplyr::select(climate_id, residuals_scaled, year)

reference_grid_shp <- shapefile("data/reference_grid.shp")

dat_grid_residuals <- st_as_sf(reference_grid_shp) %>%
  right_join(residuals, by = "climate_id") %>%
  filter(year %in% 1987:2016)

countries <- read_sf("data/countries_europe.shp")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, projection(countries))
world <- st_crop(world, st_bbox(countries) + c(-0.05, -0.05, 0.01, 0.01) * as.double(st_bbox(countries)))

p <- ggplot() +
  geom_sf(data = world, color = "black", fill = "lightgray") +
  geom_sf(data = dat_grid_residuals, 
          aes(fill = ifelse((residuals_scaled * 100) > 100, 100, ifelse((residuals_scaled * 100) < -100, -100, residuals_scaled * 100))), col = NA) +
  facet_wrap(~year) +
  scale_fill_gradient2(low = "#4885C1", mid = "white", high = "#AE3A4E", limits = c(-100, 100)) +
  theme_linedraw() +
  theme(panel.spacing = unit(0, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(3, "cm"),
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_colorbar(title = "Deviance in mortality from long-term trend (%)\n", 
                               title.position = "top",
                               title.align = 0.5)) +
  coord_sf(expand = FALSE) +
  geom_text(data = data.frame(year = 1987:2016), aes(x = Inf, y = Inf, label = year),
            hjust = 1.25, vjust = 5, size = 2.5)

ggsave("results/mortality_pulses_maps.pdf", p, width = 7.5, height = 8)


