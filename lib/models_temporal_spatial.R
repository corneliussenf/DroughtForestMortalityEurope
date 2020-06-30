
# Libraries ---------------------------------------------------------------

library(sf)
library(tidyverse)
library(raster)
library(lme4)
library(optimx)
library(parameters)
library(patchwork)
library(brms)
library(strucchange)

# Choose index and load data ----------------------------------------------

#index <- "VPD"
#index <- "PRECIP"
index <- "CWB"

climate_residuals <- read.csv(paste0("data/", index, "_residuals.csv"), stringsAsFactors = FALSE)

climate_residuals <- climate_residuals %>% mutate_at(.vars = vars(lag, month), as.integer)

### Remove grid cells with no forests

ecoregions_exclude <- c("Kola Peninsula tundra", "Scandinavian Montane Birch forest and grasslands")

climate_residuals <- climate_residuals %>%
  filter(!(ecoregions_name %in% ecoregions_exclude))

### Standardize drought index

climate_residuals <- climate_residuals %>%
  group_by(climate_id, lag, month) %>%
  mutate(index_scaled = (index - mean(index)) / sd(index)) %>%
  ungroup()

# Fit all combinations of lags/months -------------------------------------

lag_month <- expand_grid(lag = 1:6, month = 3:8)

fit_global_spatial_temporal <- vector("list", nrow(lag_month))

for (i in 1:nrow(lag_month)) {
  
  print(paste0(i, " out of ", nrow(lag_month)))
  
  fit_global_spatial_temporal[[i]] <- lmer(residuals_scaled ~ index_scaled + (1 + index_detrend_scaled | climate_id) + (1 + index_detrend_scaled | year),
                                   data = climate_residuals %>% filter(., month == lag_month[i, "month"][[1]] & lag == lag_month[i, "lag"][[1]]),
                                   control = lmerControl(optimizer = "optimx", optCtrl = list(method = "nlminb")))
  
  
}

save(fit_global_spatial_temporal, file = paste0("temp/fit_global_spatial_temporal_", index, ".RData"))
load(file = paste0("temp/fit_global_spatial_temporal_", index, ".RData"))

# Evaluate models and select lag/smoothing combination -------------------------

aic_spatial_temporal <- fit_global_spatial_temporal %>% 
  map(AIC) %>%
  map(~ data.frame(aic = .)) %>%
  set_names(paste(lag_month$lag, lag_month$month, sep = ".")) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("lag", "month"), "\\.")

write_csv(aic_spatial_temporal, paste0("results/spatial_temporal/aic_spatial_temporal_", index, ".csv"))

p_aic_spatial_temporal <- ggplot(aic_spatial_temporal, aes(x = str_pad(lag, 2, side = "left", pad = "0"),
                                                           y = str_pad(month, 2, side = "left", pad = "0"),
                                                           fill = aic)) +
  geom_tile() +
  geom_point(aes(size = ifelse(aic %in% sort(aic)[1], 1, NA)), show.legend = FALSE, shape = 4, col = "white") +
  scale_fill_viridis_c() +
  theme_linedraw() +
  theme(panel.spacing = unit(c(0, 0, 0 ,0), "cm"),
        legend.position = "right") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Temporal averaging (months)", y = "Observation month", fill = "AIC", 
       title = index) +
  guides(fill =  guide_colorbar(barheight = unit(4, "cm"), barwidth = unit(0.4, "cm"), title.position = "top"))

ggsave(paste0("results/spatial_temporal/aic_lag_month_spatial_temporal_", index, ".pdf"), p_aic_spatial_temporal, width = 3.5, height = 3)

fit_global_spatial_temporal_select <- fit_global_spatial_temporal[[which.min(aic_spatial_temporal$aic)]]

MuMIn::r.squaredGLMM(fit_global_spatial_temporal_select)

isSingular(fit_global_spatial_temporal_select)

climate_residuals_select <- climate_residuals %>% 
  filter(., month == aic_spatial_temporal[which.min(aic_spatial_temporal$aic), "month"] & 
           lag == aic_spatial_temporal[which.min(aic_spatial_temporal$aic), "lag"])

# Reestimate using Bayesian methods ---------------------------------------

aic_spatial_temporal <- read_csv(paste0("results/spatial_temporal/aic_spatial_temporal_", index, ".csv"))

climate_residuals_select <- climate_residuals %>% 
  filter(., month == aic_spatial_temporal[which.min(aic_spatial_temporal$aic), "month"][[1]] & 
           lag == aic_spatial_temporal[which.min(aic_spatial_temporal$aic), "lag"][[1]])

fit_brm_null <- brm(residuals_scaled ~ 1 + (1 | climate_id) + (1 | year),
                    data = climate_residuals_select,
                    family = "exgaussian",
                    chains = 4,
                    iter = 4000,
                    cores = 4)

fit_brm <- brm(residuals_scaled ~ index_scaled + (1 + index_scaled | climate_id) + (1 + index_scaled | year),
               data = climate_residuals_select,
               family = "exgaussian",
               chains = 4,
               iter = 4000,
               cores = 4)

fit_brm_wiggle <- brm(residuals_scaled ~ s(index_scaled) + (1 + index_scaled | climate_id) + (1 + index_scaled | year),
                      data = climate_residuals_select,
                      family = "exgaussian",
                      chains = 4,
                      iter = 4000,
                      cores = 4)


pp_check(fit_brm_wiggle, type = "dens_overlay") + xlim(-1, 3)

save(fit_brm_null, file = paste0("temp/fit_brm_spatial_temporal_null_", index, ".RData"))
save(fit_brm, file = paste0("temp/fit_brm_spatial_temporal_", index, ".RData"))
save(fit_brm_wiggle, file = paste0("temp/fit_brm_spatial_temporal_wiggle_", index, ".RData"))

load(file = paste0("temp/fit_brm_spatial_temporal_", index, ".RData"))
load(file = paste0("temp/fit_brm_spatial_temporal_null_", index, ".RData"))
load(file = paste0("temp/fit_brm_spatial_temporal_wiggle_", index, ".RData"))

pp_check(fit_brm_wiggle, type = "dens_overlay") +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1.02, 0),
        legend.justification = c(1.1, 0),
        legend.background = element_blank(),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.4), "cm"))

pp_check(fit_brm_wiggle, type = "stat", stat = "mean", nsamples = 500)
pp_check(fit_brm_wiggle, type = "stat", stat = "median", nsamples = 500)
pp_check(fit_brm_wiggle, type = "stat", stat = "min", nsamples = 500)
pp_check(fit_brm_wiggle, type = "stat", stat = "max", nsamples = 500)
pp_check(fit_brm_wiggle, type = "stat", stat = "sd", nsamples = 500)

fit_brm <- add_criterion(fit_brm, criterion = "loo", nsamples = 1000)
fit_brm_wiggle <- add_criterion(fit_brm_wiggle, criterion = "loo", nsamples = 1000)
fit_brm_null <- add_criterion(fit_brm_null, criterion = "loo", nsamples = 1000)

loo_compare(fit_brm_null, fit_brm, fit_brm_wiggle)

summary(fit_brm_wiggle)

bayes_R2(fit_brm_null, nsamples = 2000)
bayes_R2(fit_brm, nsamples = 2000)
bayes_R2(fit_brm_wiggle, nsamples = 2000)

# Calculate probability of excess mortality -------------------------------

predictions_full <- posterior_predict(fit_brm_wiggle,
                                      newdata = data.frame(index_scaled = seq(-3, 3, 0.05)),
                                      re.form = NA)

predictions_full_summary <- predictions_full %>%
  as.data.frame(.) %>%
  mutate(draws = 1:n()) %>%
  gather(., key = key, value = prediction, -draws) %>%
  mutate(CWB = rep(seq(-3, 3, 0.05), each = max(draws))) %>%
  dplyr::select(-key) %>%
  mutate(draws_group = cut(draws, seq(0, max(draws), length.out = 100))) %>%
  group_by(CWB, draws_group) %>%
  summarise(lt_0.00 = mean(prediction <= 0),
            gt_0.00 = mean(prediction > 0),
            gt_0.25 = mean(prediction > 0.25),
            gt_0.50 = mean(prediction > 0.5),
            gt_1.00 = mean(prediction > 1)) %>%
  gather(., key = excess, value = p, -CWB, -draws_group) %>%
  group_by(CWB, excess) %>%
  summarize(p_mean = mean(p),
            p_lower = quantile(p, 0.025),
            p_upper = quantile(p, 0.975)) %>%
  ungroup(.) %>%
  mutate(excess = factor(excess, 
                         levels = c("lt_0.00", "gt_0.00", "gt_0.25", "gt_0.50", "gt_1.00"), 
                         labels = c("No excess mortality", "Excess mortality", "> 25 %", "> 50 %", "> 100 %")))

### Identify threshold

predictions_full_bp <- predictions_full %>%
  as.data.frame(.) %>%
  mutate(draws = 1:n()) %>%
  gather(., key = key, value = prediction, -draws) %>%
  mutate(CWB = rep(seq(-3, 3, 0.05), each = max(draws))) %>%
  dplyr::select(-key) %>%
  mutate(draws_group = cut(draws, seq(0, max(draws), length.out = 100))) %>%
  mutate(draws_group = as.integer(draws_group)) %>%
  group_by(CWB, draws_group) %>%
  summarise(lt_000 = mean(prediction <= 0),
            gt_000 = mean(prediction > 0),
            gt_025 = mean(prediction > 0.25),
            gt_050 = mean(prediction > 0.5),
            gt_100 = mean(prediction > 1)) %>%
  ungroup() %>%
  gather(., key = excess, value = p, -CWB, -draws_group) %>%
  split(list(.$draws_group, .$excess)) %>%
  map(~ .$CWB[breakpoints(p ~ CWB, data = ., breaks = 1)$breakpoints]) %>%
  map(~ data.frame(CWB_bp = .)) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("draws_group", "excess"), "\\.") %>%
  group_by(excess) %>%
  summarize(CWB_bp_mean = mean(CWB_bp, na.rm = TRUE),
            CWB_bp_lower = quantile(CWB_bp, 0.025, na.rm = TRUE),
            CWB_bp_upper = quantile(CWB_bp, 0.975, na.rm = TRUE)) %>%
  mutate(excess = factor(excess, 
                         levels = c("lt_000", "gt_000", "gt_025", "gt_050", "gt_100"), 
                         labels = c("No excess mortality", "Excess mortality", "> 25 %", "> 50 %", "> 100 %")))

p <- ggplot(data = predictions_full_summary %>% filter(excess != "No excess mortality")) +
  geom_step(aes(x = CWB, y = p_mean, col = excess)) +
  geom_text(data = predictions_full_bp %>% filter(excess != "No excess mortality"), 
            aes(x = CWB_bp_mean + 0.05, y = c(0.7, 0.9, 0.8, 0.6), 
                label = paste0(round(CWB_bp_mean, 2), " (", round(CWB_bp_lower, 2), " - ", round(CWB_bp_upper, 2), ")")), 
            hjust = 0, col = "black", size = 2) +
  geom_segment(data = predictions_full_bp %>% filter(excess != "No excess mortality"), 
               aes(x = CWB_bp_mean, xend = CWB_bp_mean, 
                   y = 0, yend = c(0.7, 0.9, 0.8, 0.6),
                   col = excess), linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1.02, 1.05),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.4), "cm")) +
  labs(x = "Climatic water balance (z-scores)", y = "Probability", col = NULL, fill = NULL) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, name = "RdBu")[c(4, 3, 2, 1)]) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, name = "RdBu")[c(4, 3, 2, 1)])

ggsave(paste0("results/spatial_temporal/thresholds_", index, ".pdf"), p, width = 3.5, height = 3.5)

write_csv(predictions_full_bp, paste0("results/spatial_temporal/thresholds_", index, ".csv"))

### Final plot

p <- ggplot(data = predictions_full_summary) +
  geom_step(aes(x = CWB, y = p_mean, col = excess)) +
  pammtools::geom_stepribbon(aes(x = CWB, ymin = p_lower, ymax = p_upper, fill = excess), alpha = 0.25) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1.02, 1.05),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.4), "cm")) +
  labs(x = "Climatic water balance (z-scores)", y = "Probability", col = NULL, fill = NULL) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, name = "RdBu")[c(8, 4, 3, 2, 1)]) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, name = "RdBu")[c(8, 4, 3, 2, 1)])

ggsave(paste0("results/spatial_temporal/probability_excess_mortality_", index, ".pdf"), p, width = 3, height = 3)

# Numbers

predictions_full_summary %>%
  filter(CWB %in% c(-1, -2, -3) & excess == "Excess mortality")

# Smoothing term ----------------------------------------------------------

predictions_expectation <- pp_expect(fit_brm_wiggle,
                                     newdata = data.frame(index_scaled = seq(-3, 3, 0.1)),
                                     re.form = NA)

predictions_expectation_summary <- data.frame(CWB = seq(-3, 3, 0.1),
                                              pred_mean = apply(predictions_expectation, 2, mean),
                                              pred_sd = apply(predictions_expectation, 2, sd),
                                              pred_median = apply(predictions_expectation, 2, quantile, 0.5),
                                              pred_lower = apply(predictions_expectation, 2, quantile, 0.75),
                                              pred_upper = apply(predictions_expectation, 2, quantile, 0.25),
                                              pred_llower = apply(predictions_expectation, 2, quantile, 0.025),
                                              pred_uupper = apply(predictions_expectation, 2, quantile, 0.975))

predictions_expectation <- predictions_expectation %>%
  as.data.frame(.) %>%
  mutate(draw = 1:n()) %>%
  gather(key = key, value = pred, -draw) %>%
  mutate(key = rep(seq(-3, 3, 0.1), each = max(draw)))

# Summarized data points before plotting

climate_residuals_select_summary <- climate_residuals_select %>%
  group_by(cwb_cut = cut(index_detrend_scaled, 
                         breaks = c(-4, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 4), 
                         labels = c(-2.9, -2, -1, 0, 1, 2, 2.9))) %>%
  summarise(mn = mean(residuals_scaled),
            sd = sd(residuals_scaled),
            min = quantile(residuals_scaled, 0.01),
            max = quantile(residuals_scaled, 0.99)) %>%
  mutate(CWB = as.double(as.character(cwb_cut)))

p1 <- ggplot() +
  geom_errorbar(data = climate_residuals_select_summary,
                aes(x = CWB, ymin = mn - sd,  ymax = mn + sd), 
                col = "grey", width = 0) +
  geom_point(data = climate_residuals_select_summary,
             aes(x = CWB, y = mn), 
             col = "grey") +
  geom_line(data = predictions_expectation_summary,
            aes(x = CWB, y = pred_mean), col = "#AE3A4E") +
  geom_ribbon(data = predictions_expectation_summary,
              aes(x = CWB, ymin = pred_lower, ymax = pred_upper), 
              fill = "#AE3A4E", alpha = 0.25) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 2), label = function(x) x * 100) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_linedraw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.4), "cm")) +
  labs(x = "Climatic water balance (z-scores)", y = "Deviation in canopy mortality\nfrom long-term trend (%)")

ggsave(paste0("results/spatial_temporal/response_curve_mortality_", index, ".pdf"), p1, width = 3, height = 2.75)

# Full data points

p2 <- ggplot() +
  geom_point(data = climate_residuals_select %>%
               sample_frac(0.1),
             aes(x = index_detrend_scaled, y = residuals_scaled), 
             col = "grey", alpha = 0.1) +
  geom_line(data = predictions_expectation_summary,
            aes(x = CWB, y = pred_mean), col = "#AE3A4E") +
  geom_ribbon(data = predictions_expectation_summary,
              aes(x = CWB, ymin = pred_llower, ymax = pred_uupper), 
              fill = "#AE3A4E", alpha = 0.25) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 3), label = function(x) x * 100) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_linedraw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.4), "cm")) +
  labs(x = "Climatic water balance (z-scores)", y = "Deviation in canopy mortality\nfrom long-term trend (%)")

ggsave(paste0("results/spatial_temporal/response_curve_mortality_", index, "_full.pdf"), p2, width = 3, height = 2.75)

# Both combined

p <- p2 + p1 + plot_layout(ncol = 2)

ggsave("results/spatial_temporal/response_curve_mortality_CWB_combined.pdf", plot = p, width = 7.5, height = 3.5)

# Numbers

pred <- pp_expect(fit_brm_wiggle,
          newdata = data.frame(index_scaled = seq(-3, 3, 1)),
          re.form = NA)

pred <- data.frame(CWB = seq(-3, 3, 1),
                   mean = round(apply(pred, 2, mean) * 100, 2),
                   lower = round(apply(pred, 2, quantile, 0.025) * 100, 2),
                   upper = round(apply(pred, 2, quantile, 0.975) * 100, 2))

write_csv(pred, paste0("results/average_change_", index, ".csv"))

# Get model parameters ----------------------------------------------------

#load(file = paste0("temp/fit_brm_", index, ".RData"))

posterior <- as.data.frame(fit_brm_wiggle)

fixedeffect <- posterior %>%
  summarize(mean = mean(bs_sindex_scaled_1),
            lower = quantile(bs_sindex_scaled_1, 0.025),
            upper = quantile(bs_sindex_scaled_1, 0.975)) %>%
  mutate(index = index)

write_csv(fixedeffect, paste0("results/spatial_temporal/fixedeffect_", index, ".csv"))

# Plot spatial variability in effect size ---------------------------------------

climate_ids_names <- grep(glob2rx("r_climate_id[*,index_scaled]"), names(posterior), value = TRUE)
climate_ids <- climate_ids_names %>% sub("r_climate_id\\[", "", .) %>% sub(",index_scaled\\]", "", .) %>% as.integer(.)

ranefs_climate_id <- data.frame(climate_id = climate_ids,
                                effect = apply(posterior[, climate_ids_names], 2, mean))
rownames(ranefs_climate_id) <- NULL

write_csv(ranefs_climate_id, paste0("results/spatial_temporal/randefs_cliamteid_", index, ".csv"))

reference_grid_shp <- read_sf("data/reference_grid.shp")

countries <- read_sf("data/countries_europe.shp")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, projection(countries))
st_crs(world) <- st_crs(countries)
world <- st_crop(world, st_bbox(countries) + c(-0.05, -0.05, 0.01, 0.01) * as.double(st_bbox(countries)))

ranefs_climate_id <- reference_grid_shp %>%
  st_as_sf() %>%
  right_join(ranefs_climate_id %>% mutate(climate_id = as.integer(climate_id)), 
             by = "climate_id")

maxlegend <- c(-1 * max(abs(ranefs_climate_id$effect)), max(abs(ranefs_climate_id$effect)))

p_effect_cliamte_id <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = ranefs_climate_id, aes(fill = effect), col = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  theme_linedraw() +
  #scale_fill_gradient2(low = "#33BBEE", mid = "white", high = "#CC3311", limits = maxlegend) +
  scale_fill_gradient2(low = "#4885C1", mid = "white", high = "#AE3A4E", limits = maxlegend) +
  labs(fill = bquote("Spatial variation in slope ("*beta[i]*")")) +
  theme(panel.spacing = unit(0, "cm"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(hjust = 1),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.title.align = 0.5,
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  coord_sf(expand = FALSE) +
  guides(fill = guide_colorbar(title.position = "top"))

if (index != "VPD") {
  p_effect_cliamte_id <- p_effect_cliamte_id +
    scale_fill_gradient2(low = "#AE3A4E", mid = "white", high = "#4885C1", limits = maxlegend)
}

ggsave(paste0("results/spatial_temporal/ranef_climate_id_", index, ".pdf"), p_effect_cliamte_id, width = 4.5, height = 5.5)

# Plot temporal variability in effect size ----------------------------------------

years_names <- grep(glob2rx("r_year[*,index_scaled]"), names(posterior), value = TRUE)
years <- years_names %>% sub("r_year\\[", "", .) %>% sub(",index_scaled\\]", "", .) %>% as.integer(.)

ranefs_years <- data.frame(year = years,
                           effect = apply(posterior[, years_names], 2, median),
                           lower = apply(posterior[, years_names], 2, quantile, 0.25),
                           upper = apply(posterior[, years_names], 2, quantile, 0.75),
                           llower = apply(posterior[, years_names], 2, quantile, 0.1),
                           uupper = apply(posterior[, years_names], 2, quantile, 0.9))
rownames(ranefs_years) <- NULL

write_csv(ranefs_years, paste0("results/spatial_temporal/randefs_years_", index, ".csv"))

p_effect_years <- ggplot(data = ranefs_years) +
  geom_errorbar(aes(x = year, ymin = lower, ymax = upper), width = 0) +
  geom_point(aes(x = year, y = effect, col = effect), size = 2) +
  theme_linedraw() +
  geom_hline(yintercept = 0) +
  scale_color_gradient2(low = "#AE3A4E", mid = "white", high = "#4885C1", limits = c(-0.12, 0.12)) +
  labs(x = "Year", y = bquote("Temporal variation in slope ("*beta[t]*")"), col = bquote("Temporal variation in slope ("*beta[t]*")")) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.title.align = 0.5,
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        plot.margin = unit(c(0.1, 0.1, 0.5, 0.5), "cm")) +
  guides(color = guide_colorbar(title.position = "top"))

ggsave(paste0("results/spatial_temporal/ranef_years_", index, ".pdf"), p_effect_years, width = 5.5, height = 5.5)

