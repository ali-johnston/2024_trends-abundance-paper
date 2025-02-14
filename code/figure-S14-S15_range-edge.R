library(tidyverse)
library(lme4)
library(arrow)
library(fs)
library(ggplot2)

data_dir <- "data"
outputs_dir <- "outputs"
figures_dir <- path("figures", "figure-S14-S15")
dir_create(figures_dir)


# load data ----

# trends estimates
trends <- path(data_dir, "ebird-trends_2007-2021.parquet") |>
  read_parquet() |>
  mutate(log_abd = log10(abd),
         log_distance_to_edge_km = log10(distance_to_edge_km),
         breeding_biome = factor(breeding_biome),
         species_code = factor(species_code),
         weight = abd_ppy_sd / mean(abd_ppy_sd))

# minimum and maximum distance to edge for each species and each biome
# used for predictions and plots
dist_range <- trends |>
  dplyr::select(species_code, breeding_biome, distance_to_edge_km) |>
  group_by(species_code, breeding_biome) |>
  summarise(max_dist = max(distance_to_edge_km),
            min_dist = min(distance_to_edge_km),
            .groups = "drop") |>
  mutate(min_log_dist = log10(min_dist), max_log_dist = log10(max_dist))
dist_range_biome <- dist_range |>
  dplyr::select(breeding_biome, min_log_dist, max_log_dist) |>
  group_by(breeding_biome) |>
  summarise(min_log_dist = min(min_log_dist),
            max_log_dist = max(max_log_dist), .groups = "drop")

# Pick up model results from bam
mod_tag <- "log10_distance_to_edge_km_40_run2"
results_loc1 <- path(outputs_dir, paste0("species_coefs_bam_linear_", mod_tag, ".csv"))
results_loc2 <- path(outputs_dir, paste0("species_coefs_lm_linear_", mod_tag, ".csv"))

spec_coef <- read_csv(results_loc1) |>
  bind_rows(read_csv(results_loc2)) |>
  mutate(sig = ifelse(p_val < 0.05, "s", "ns")) |>
  mutate(sig_fac = factor(sig, levels = c("s", "ns"), ordered = TRUE)) |>
  rename(species_code = species_code...5) |>
  dplyr::select(-species_code...10) |>
  merge(dist_range, by = "species_code") |>
  mutate(effect_size = est*(max_log_dist - min_log_dist)) |>
  mutate(intercept0 = int_est + est*min_log_dist) |>
  mutate(intercept0.5 = int_est + est*mean(c(min_log_dist, max_log_dist)))


# mean effect sizes per biome
biome_coef <- spec_coef |>
  dplyr::select(breeding_biome, intercept0, effect_size) |>
  rename(est = effect_size, int_est = intercept0) |>
  group_by(breeding_biome) |>
  summarise_all(mean) |>
  mutate(min_log_dist = 0, max_log_dist = 1)

# figure S13: histogram of slopes ----

# max slope for limits
max_slope <- ceiling(max(abs(spec_coef$est))*10)/10
max_slope <- 10

ggplot(spec_coef, aes(x = est, fill = sig_fac)) +
  geom_histogram(breaks = seq(-1*max_slope, max_slope, length.out = 41),
                 colour = "white", show.legend = FALSE) +
  scale_fill_manual(values = c(alpha("steelblue3", 0.8), alpha("steelblue3", 0.4))) +
  theme_bw() +
  geom_vline(xintercept = 0) +
  ylab("Number of species") +
  xlab("Linear model slope") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(path(figures_dir, "figure-S14_linear-model_slope-histogram.png"),
       width = 14, height = 10, units = "cm")
ggsave(path(figures_dir, "figure-S14_linear-model_slope-histogram.tif"),
       width = 14, height = 10, units = "cm")


# figure S15: trend vs. distance to edge relationships from random effects ----

# summarise significance and direction
spec_coef <- spec_coef |> 
  mutate(direction = ifelse(est>0, "pos", "neg"))
spec_coef |>
  dplyr::select(sig_fac, direction) |>
  table()
#        direction
# sig_fac neg pos
#      s  195 121
#      ns  95  84

sum(spec_coef$direction == "neg") / nrow(spec_coef)
# [1] 0.5858586
sum(spec_coef$direction == "neg" & spec_coef$sig_fac == "s") / nrow(spec_coef)
# [1] 0.3939394

nu <- function(x){as.numeric(as.character(x))}
plot_1_line <- function(data, col = alpha("black", 0.1), ...){
  x <- seq(0, 1, by = 0.1)
  x2 <- seq(data$min_log_dist[1], data$max_log_dist[1], length.out = 11)
  lines(x, x2*data$est + data$int_est, col = col, ...)
}

plot_multi <- function(data, col = alpha("black", 0.1)){
  for(i in 1:nrow(data)){
    plot_1_line(data[i,], col = col)
  }
}

# change breeding biome names to match paper
biomes <- names(table(spec_coef$breeding_biome))
biome_names_tidy <- c("Arctic tundra", "Aridland", "Forest",
                      "Grassland", "Habitat generalists", "Wetland & Coast")
tiff(path(figures_dir, "figure-S15_linear-model_slopes.tif"),
     width = 14, height = 10, units = "cm", pointsize = 9, res = 600)

par(mfrow=c(2, 3), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,5,2,0))
biomes <- names(table(spec_coef$breeding_biome))

for (i in seq_along(biomes)) {
  ylabel <- ifelse(i %in% c(1, 4), "Population trend", "")
  plot(0, 0, col = "white", ylim = c(-12, 12), xlim = c(0, 1),
       main = "", xaxt = "n", xlab = "", yaxt = "n", ylab = ylabel,
       xpd = NA, cex.lab = 1.5)
  if(i %in% c(1,4)) axis(side = 2, at = c(-10, -5, 0, 5, 10), cex.axis = 1.5)
  if(i > 3) axis(side = 1, at = c(0, 1), labels = c("0", "max"), cex.axis = 1.5)
  abline(h = 0, col = alpha("firebrick", 0.42), lwd = 2)
  
  # lines for each species in the biome
  plot_multi(filter(spec_coef, breeding_biome == biomes[i]))
  
  if(i == 5) text(x = 0.5, y = -17,
                  labels = "Relative distance to edge within species", xpd = NA, cex = 1.5)
  text(x = 0, y = 11, pos = 4, labels = biome_names_tidy[i], cex = 1.5)
  
  # add average biome line
  biome_avg <- filter(biome_coef, breeding_biome == biomes[i])
  plot_1_line(biome_avg, col = "black", lwd = 2)
}

dev.off()

# biome summaries

biome_coef |> arrange(est)
#   breeding_biome     int_est      est min_log_dist max_log_dist
#   <chr>                <dbl>    <dbl>        <dbl>        <dbl>
# 1 Aridlands            0.627 -0.429              0            1
# 2 Forest              -0.346 -0.161              0            1
# 3 Wetlands & Coasts   -1.39  -0.0667             0            1
# 4 Habitat Generalist  -0.555 -0.0616             0            1
# 5 Grassland           -1.75  -0.00166            0            1
# 6 Arctic Tundra       -2.61   0.0606             0            1


spec_coef |> 
  mutate(direction = ifelse(est < 0, "neg", "pos")) |>
  mutate(sig_neg = ifelse(direction == "neg" & sig == "s", 1, 0)) |>
  mutate(neg_bin = ifelse(est < 0, 1, 0)) |>
  dplyr::select(breeding_biome, sig_neg, neg_bin) |>
  group_by(breeding_biome) |>
  summarise_all(mean) |>
  arrange(neg_bin)

#   breeding_biome     sig_neg neg_bin
#   <chr>                <dbl>   <dbl>
# 1 Habitat Generalist   0.415   0.463
# 2 Arctic Tundra        0.167   0.467
# 3 Wetlands & Coasts    0.303   0.582
# 4 Forest               0.433   0.596
# 5 Grassland            0.478   0.609
# 6 Aridlands            0.487   0.671