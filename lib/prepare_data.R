
# Libraries ---------------------------------------------------------------

setwd("LandTrendr/drought/")

library(sf)
library(tidyverse)
library(raster)
library(Rcpp)
library(fasterize)
library(lme4)
library(optimx)
library(parameters)
library(patchwork)

sourceCpp("lib/fasttable.cpp")

# Get refernce grid and country outlines -----------------------------------

reference_grid <- raster("data/climatevars/CWB/1986.01.tif")
reference_grid[!is.na(reference_grid)] <- 1:cellStats(!is.na(reference_grid), sum)
reference_grid_shp <- rasterToPolygons(reference_grid)
reference_grid_shp <- spTransform(reference_grid_shp, CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
names(reference_grid_shp) <- "climate_id"

rgdal::writeOGR(reference_grid_shp, "data", "reference_grid", driver = "ESRI Shapefile")

writeRaster(reference_grid, "data/climate_id_reference_grid.tif", datatype = "INT2U", overwrite = TRUE)

countries <- read_sf("data/countries_europe.shp")

country_master <- read_csv("data/countries_master.csv") %>%
  filter(country_name_short != "russia")

# Aggregate disturbances to climate grid ----------------------------------

for (i in 1:nrow(country_master)) {
  
  cntr <- country_master[i, "country_name_short"][[1]]
  
  print(cntr)
  
  disturbance <- raster(paste0("../prediction/", cntr, "/disturbance_year_filtered_", cntr, ".tif"))
  forest <- raster(paste0("../prediction/", cntr, "/prediction_forestcover_", cntr, ".tif"))
  
  ext <- as(extent(disturbance), 'SpatialPolygons')
  proj4string(ext) <- projection(disturbance)
  
  grid_sel <- st_intersection(st_as_sf(reference_grid_shp), st_as_sf(ext))
  grid_sel_ras <- fasterize(grid_sel, disturbance, field = "climate_id")
  grid_values <- values(grid_sel_ras)
  
  dat <- data.frame(climate_id = grid_values,
                    disturbance = values(disturbance),
                    country = cntr) %>%
    na.omit(.) %>%
    group_by(climate_id, year = disturbance) %>%
    summarize(disturbance_ha = n() * 0.09,
              country = unique(country)) %>%
    ungroup(.)
  
  forest <- data.frame(climate_id = grid_values,
                       forest = values(forest)) %>%
    filter(!is.na(forest)) %>%
    group_by(climate_id) %>%
    summarize(forest_ha = sum(forest == 1, na.rm = TRUE) * 0.09,
              land_ha = n() * 0.09) %>%
    ungroup(.)
  
  dat <- dat %>% 
    left_join(forest, by = "climate_id")
  
  write_csv(dat, paste0("data/disturbances/grid_aggregate_", cntr, ".csv"))
  
}

# Ecoregions and biomes ----------------------------------------------------

ecoregions <- read_sf("data/ecoregions/terrestrial_ecoregions_olson.shp")

ecoregions_grid <- fasterize(ecoregions %>% st_transform(., projection(reference_grid)), 
                             reference_grid, field = "ECO_ID")

ecoregions_dat <- data.frame(climate_id = values(reference_grid),
                             ecoregions = values(ecoregions_grid)) %>%
  na.omit(.) 

ecoregions_dat <- ecoregions_dat %>%
  left_join(ecoregions %>% 
              as_tibble(.) %>% 
              group_by(ecoregions = ECO_ID) %>% 
              summarize(ecoregions_name = unique(ECO_NAME),
                        biome = unique(BIOME)), by = "ecoregions")

# Load aggregated data ----------------------------------------------------

reference_grid <- raster("data/climatevars/cliamte_id_reference_grid.tif")
reference_grid_shp <- rasterToPolygons(reference_grid)
reference_grid_shp <- spTransform(reference_grid_shp, CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
names(reference_grid_shp) <- "climate_id"

countries <- read_sf("data/countries_europe.shp")

country_master <- read_csv("data/countries_master.csv") %>%
  filter(country_name_short != "russia")

names <- list.files("data/disturbances", pattern = glob2rx("*grid_aggregate*.csv")) %>%
  strsplit(., "_") %>%
  map(., ~ substring(.[4], 1, nchar(.[4]) - 4)) %>%
  unlist()

dat <- list.files("data/disturbances", pattern = glob2rx("*grid_aggregate*.csv"), full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()

dat_forest <- dat %>%
  group_by(climate_id, country) %>%
  summarize(forest_ha = unique(forest_ha),
            land_ha = unique(land_ha)) %>%
  group_by(climate_id) %>%
  summarise(forest_ha = sum(forest_ha),
            land_ha = sum(land_ha)) %>%
  ungroup() 

dat <- dat %>%
  group_by(climate_id, year) %>%
  summarize(disturbance_ha = sum(disturbance_ha, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(dat_forest, by = "climate_id") %>%
  mutate(rate = disturbance_ha / forest_ha)

dat <- dat %>%
  split(.$climate_id) %>%
  map(~ right_join(., data.frame(climate_id = unique(.$climate_id),
                                 forest_ha = unique(.$forest_ha),
                                 land_ha = unique(.$land_ha),
                                 year = 1986:2018), 
                   by = c("climate_id", "year", "forest_ha", "land_ha"))) %>%
  map(~ mutate_at(., .vars = vars(disturbance_ha, rate), .funs = function(x) ifelse(is.na(x), 0, x))) %>%
  bind_rows()

write_csv(dat, "temp/dat.csv")

dat_grid <- st_as_sf(reference_grid_shp) %>%
  right_join(dat, by = "climate_id") %>%
  left_join(ecoregions_dat, by = "climate_id")

ggplot() +
  geom_sf(data = dat_grid, aes(fill = rate), col = NA) +
  facet_wrap(~year) +
  scale_fill_viridis_c(limits = c(0, 0.03))

ggplot() +
  geom_sf(data = dat_grid %>% filter(year == 2000), aes(fill = ecoregions_name), col = NA)

# Filter for years

dat <- dat %>% filter(year %in% 1987:2016)

# Calculate residuals (i.e., disturbance pulses) ---------------------------

# Example figure

climate_id_exp <- sample(unique(dat$climate_id), 1)
# Good one: 25296
# climate_id_exp <- 25296

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
  labs(y = "Canopy mortality (% per year)", col = NULL) +
  scale_color_manual(values = "lightgrey")

p_exp2 <- ggplot(dat_exp, aes(x = year, y = residuals_relative * 100, fill = residuals_relative)) +
  geom_bar(stat = "identity") +
  theme_linedraw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.margin = unit(c(0, 0.1, 0.5, 0.5), "cm")) +
  labs(y = "Percent change\nfrom trend", x = "Year") +
  scale_fill_gradient2(low = "#33BBEE", mid = "white", high = "#CC3311") +
  geom_hline(yintercept = 0)
  
p_exp <- p_exp1 + p_exp2 + plot_layout(ncol = 1, height = c(0.6, 0.4)) 

p_exp

ggsave(paste0("results/residuals_example_", climate_id_exp, ".pdf"), p_exp, width = 3.5, height = 3.5)

# Calculate for each cell

residuals <- dat %>%
  split(.$climate_id) %>%
  map(~ mutate(., disturbance_n = disturbance_ha / 0.09,
               forest_n = forest_ha / 0.09)) %>%
  map(~ glm(cbind(disturbance_n, forest_n) ~ year, data = ., family = binomial(link = "logit"))) %>%
  map(~ data.frame(year = .$data$year,
                   rate = .$data$rate,
                   fitted = fitted(.x))) %>%
  map(~ mutate(.,
               residuals = rate - fitted,
               residuals_scaled = residuals / fitted)) %>%
  bind_rows(.id = "climate_id") %>%
  mutate(climate_id = as.integer(climate_id)) %>%
  left_join(ecoregions_dat, by = "climate_id") %>%
  filter(!is.na(ecoregions))

# Temporal

p <- ggplot(data = residuals %>%
         group_by(year) %>%
         summarize(residuals_scaled_mn = mean(residuals_scaled, na.rm = TRUE),
                   residuals_scaled_sd = sd(residuals_scaled, na.rm = TRUE),
                   residuals_scaled_se = sd(residuals_scaled, na.rm = TRUE) / sqrt(sum(!is.na(residuals_scaled))))) +
  geom_ribbon(aes(x = year, 
                  ymin = residuals_scaled_mn - residuals_scaled_sd, 
                  ymax = residuals_scaled_mn + residuals_scaled_sd), alpha = 0.5) +
  geom_line(aes(x = year, y = residuals_scaled_mn)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Disturbances pulses (standardized difference to long term trend)", col = "Ecoregions", fill = "Ecoregions") +
  scale_x_continuous(breaks = c(1995, 2010)) +
  theme_linedraw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.margin = unit(c(0, 0.1, 0.5, 0.5), "cm")) +
  labs(y = "Standardized deviation from trend", x = "Year")

ggsave("results/average_residuals_europe.pdf", p, width = 3.5, height = 3)

p <- ggplot(data = residuals %>%
         mutate(ecoregions_abbr = abbreviate(ecoregions_name)) %>%
         group_by(year, ecoregions_name, ecoregions_abbr) %>%
         summarize(residuals_scaled_mn = mean(residuals_scaled, na.rm = TRUE),
                   residuals_scaled_sd = sd(residuals_scaled, na.rm = TRUE),
                   residuals_scaled_se = sd(residuals_scaled, na.rm = TRUE) / sqrt(sum(!is.na(residuals_scaled))))) +
  geom_ribbon(aes(x = year, ymin = residuals_scaled_mn - residuals_scaled_sd, ymax = residuals_scaled_mn + residuals_scaled_sd, fill = ecoregions_name), alpha = 0.5) +
  geom_line(aes(x = year, y = residuals_scaled_mn, col = ecoregions_name)) +
  facet_wrap(~ecoregions_abbr) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  guides(colour = guide_legend(ncol = 1, keywidth = unit(0.2, "cm"), keyheight = unit(0.2, "cm")),
         fill = guide_legend(ncol = 1, keywidth = unit(0.2, "cm"), keyheight = unit(0.2, "cm"))) +
  labs(x = "Year", y = "Disturbances pulses (standardized difference to long term trend)", col = "Ecoregions", fill = "Ecoregions") +
  scale_x_continuous(breaks = c(1995, 2010)) +
  theme_linedraw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.margin = unit(c(0, 0.1, 0.5, 0.5), "cm")) +
  labs(y = "Standardized deviation from trend", x = "Year")

ggsave("results/average_residuals_europe_ecoregions.pdf", p, width = 7.5, height = 7.5)

# Spatial

dat_grid_residuals <- st_as_sf(reference_grid_shp) %>%
  right_join(residuals, by = "climate_id") %>%
  filter(year %in% 1987:2016)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, projection(countries))
world <- st_crop(world, st_bbox(countries) + c(-0.05, -0.05, 0.01, 0.01) * as.double(st_bbox(countries)))

p <- ggplot() +
  geom_sf(data = world, color = "black", fill = "lightgray") +
  geom_sf(data = dat_grid_residuals, 
          aes(fill = ifelse(residuals_scaled < 0, 0, ifelse(residuals_scaled > 3, 3, residuals_scaled))), col = NA) +
  facet_wrap(~year) +
  scale_fill_gradient(low = "white", high = "#CC3311", breaks = c(0, 1, 2, 3), labels = c("0", "100", "200", ">300")) + 
  theme_linedraw() +
  theme(panel.spacing = unit(0, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_colorbar(title = "Disturbances pulses (percent increase in mortality compared to the long-term trend)", title.position = "top")) +
  coord_sf(expand = FALSE)

ggsave("results/disturbance_anomaly_maps.pdf", p, width = 7.5, height = 8)

# Combine with drought index ----------------------------------------------

### Select index and load data

#index <- "VPD"
#index <- "CWB"
index <- "PRECIP"

climate_var <- list.files(paste0("data/climatevars/", index), pattern = glob2rx("*tif"), full.names = TRUE) %>%
  stack(.)

### Combine with residuals and integrate -> detrend -> standardize

climate_var_sel <- vector("list", length(1984:2016))

k <- 0

for (i in seq(1, nlayers(climate_var) - 11, 12)) {
  
  k <- k + 1
  
  climate_var_sel[[k]] <- subset(climate_var, i:(i+11)) %>%
    as.data.frame(.) %>%
    mutate(climate_id = values(reference_grid)) %>%
    gather(key = year, value = index, -climate_id) %>%
    separate("year", c("year", "month"), "\\.") %>%
    mutate(year = as.integer(substr(year, 2, 5))) %>%
    na.omit(.) %>%
    left_join(residuals, by = c("year", "climate_id"))
  
}

climate_residuals <- climate_var_sel %>%
  map(bind_rows) %>%
  bind_rows()

# Filter cells not located in Europe

climate_id_eur <- unique(residuals$climate_id)

climate_residuals <- climate_residuals %>%
  filter(., climate_id %in% climate_id_eur)

# Calculate integrated climate variables

climate_residuals <- climate_residuals %>%
  split(.$climate_id) %>%
  map(~ arrange(., year, month)) %>%
  map(~ mutate(.,
               index_ma01 = zoo::rollmean(index, k = 1, na.pad = TRUE, align = "right"),
               index_ma02 = zoo::rollmean(index, k = 2, na.pad = TRUE, align = "right"),
               index_ma03 = zoo::rollmean(index, k = 3, na.pad = TRUE, align = "right"),
               index_ma04 = zoo::rollmean(index, k = 4, na.pad = TRUE, align = "right"),
               index_ma05 = zoo::rollmean(index, k = 5, na.pad = TRUE, align = "right"),
               index_ma06 = zoo::rollmean(index, k = 6, na.pad = TRUE, align = "right"),
               index_ma07 = zoo::rollmean(index, k = 7, na.pad = TRUE, align = "right"),
               index_ma08 = zoo::rollmean(index, k = 8, na.pad = TRUE, align = "right"),
               index_ma09 = zoo::rollmean(index, k = 9, na.pad = TRUE, align = "right"),
               index_ma10 = zoo::rollmean(index, k = 10, na.pad = TRUE, align = "right"),
               index_ma11 = zoo::rollmean(index, k = 11, na.pad = TRUE, align = "right"),
               index_ma12 = zoo::rollmean(index, k = 12, na.pad = TRUE, align = "right"))) %>%
  bind_rows(.) %>%
  gather(key = lag, value = index, -(climate_id:biome)) %>%
  mutate(metric = substring(lag, 7, 8)) %>%
  mutate(lag = substring(lag, 9, 10))
  
# Detrend

climate_residuals <- climate_residuals %>%
  split(list(.$climate_id, .$month, .$lag)) %>%
  map(~ lm(index ~ year, data = .)) %>%
  map(~ data.frame(year = .$model$year, index_detrend = .$residuals)) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("climate_id", "month", "lag"), "\\.") %>%
  mutate(climate_id = as.integer(climate_id)) %>%
  right_join(climate_residuals, by = c("year", "month", "climate_id", "lag"))

# Standardize

climate_residuals <- climate_residuals %>%
  group_by(climate_id, month, lag) %>%
  mutate(index_detrend_scaled = index_detrend / ifelse(sd(index_detrend, na.rm = TRUE) == 0, 0, sd(index_detrend, na.rm = TRUE))) %>%
  ungroup()

# Filter to observation period

climate_residuals <- climate_residuals %>%
  filter(year %in% 1987:2016)

# Write to disc

write_csv(climate_residuals, paste0("data/", index, "_residuals.csv"))

### Plots

# Spatial

ggplot(climate_residuals %>% filter(month == "09" & lag == "09")) +
  geom_point(aes(x = index_detrend_scaled, y = residuals_scaled), alpha = 0.01, col = "red") +
  geom_smooth(aes(x = index_detrend_scaled, y = residuals_scaled), method = "lm", se = TRUE) +
  facet_wrap(~year) +
  ylim(-1, 2.5)

ggplot(climate_residuals %>% filter(month == "08" & year == "2012" & lag == "04"),
       aes(x = index_detrend_scaled, y = residuals_scaled+1)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", col = "red") +
  geom_smooth(method = "lm", col = "orange", formula = y ~ x + I(x^2)) +
  geom_smooth(method = "gam", col = "blue", formula = y ~ s(x))

summary(lm(residuals_scaled ~ index_detrend_scaled, 
           data = climate_residuals %>% 
             filter(month == "08" & year == "2012" & lag == "04")))

summary(lm(residuals_scaled ~ index_detrend_scaled + I(index_detrend_scaled^2), 
           data = climate_residuals %>% 
             filter(month == "08" & year == "2012" & lag == "04")))

# Temporal

ggplot(climate_residuals %>% filter(month == "08" & lag == "04") %>% sample_frac(., 0.25)) +
  geom_point(aes(x = index_detrend_scaled, y = residuals_scaled), alpha = 0.15) +
  facet_wrap(~ecoregions_name, scales = "free") +
  geom_smooth(aes(x = index_detrend_scaled, y = residuals_scaled), method = "lm", se = TRUE, col = "red")

