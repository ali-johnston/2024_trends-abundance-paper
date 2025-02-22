library(tidyverse)
library(lme4)
library(arrow)
library(fs)

data_dir <- "data"
outputs_dir <- "outputs"
figures_dir <- path("figures", "figure-03")
dir_create(figures_dir)


# load data ----

# trends estimates
trends <- path(data_dir, "ebird-trends_2007-2021.parquet") |>
  read_parquet() |>
  mutate(log_abd = log10(abd),
         breeding_biome = factor(breeding_biome),
         species_code = factor(species_code),
         weight = 1 / abd_ppy_var)

# minimum and maximum abundance for each species and each biome
abd_range <- trends |>
  dplyr::select(species_code, breeding_biome, abd) |>
  group_by(species_code, breeding_biome) |>
  summarise(max_abd = max(abd), min_abd = min(abd), .groups = "drop") |>
  mutate(min_log_abd = log10(min_abd), max_log_abd = log10(max_abd))

abd_range_biome <- abd_range |>
  dplyr::select(breeding_biome, min_log_abd, max_log_abd) |>
  group_by(breeding_biome) |>
  summarise(min_log_abd = min(min_log_abd),
            max_log_abd = max(max_log_abd), .groups = "drop")

# Pick up model results from bam
mod_tag <- "log10_abd_40_run2"
results_loc1 <- path(outputs_dir, paste0("species_coefs_bam_linear_", mod_tag, ".csv"))
results_loc2 <- path(outputs_dir, paste0("species_coefs_lm_linear_", mod_tag, ".csv"))

ppy_to_multiyr <- function(ppy, noyrs = 10){
  ppy_multiyr <- 100*(((100 + ppy)/100)^noyrs) - 100
  return(ppy_multiyr)
}


spec_coef <- read_csv(results_loc1) |>
  bind_rows(read_csv(results_loc2)) |>
  mutate(sig = ifelse(p_val < 0.05, "s", "ns")) |>
  mutate(sig_fac = factor(sig, levels = c("s", "ns"), ordered = TRUE)) |>
  rename(species_code = species_code...5) |>
  dplyr::select(-species_code...10) |>
  merge(abd_range, by = "species_code") |>
  mutate(effect_size = est*(max_log_abd - min_log_abd)) |>
  mutate(intercept0 = int_est + est*min_log_abd) |>
  mutate(intercept0.5 = int_est + est*mean(c(min_log_abd, max_log_abd))) |>
  mutate(ppy_min_log_abd = int_est + est*min_log_abd) |>
  mutate(ppy_max_log_abd = int_est + est*max_log_abd) |>
  mutate(multiyr_min_log_abd = ppy_to_multiyr(ppy_min_log_abd)) |>
  mutate(multiyr_max_log_abd = ppy_to_multiyr(ppy_max_log_abd)) |>
  mutate(effect_size_mutliyr = multiyr_max_log_abd - multiyr_min_log_abd)


## mean effect sizes per biome
biome_coef <- spec_coef |>
  dplyr::select(breeding_biome, intercept0, effect_size) |>
  rename(est = effect_size, int_est = intercept0) |>
  group_by(breeding_biome) |>
  summarise_all(mean) |>
  mutate(min_log_abd = 0, max_log_abd = 1)


# figure 3a: histogram of slopes ----

# max slope for limits
max_slope <- ceiling(max(abs(spec_coef$est))*10)/10
ggplot(spec_coef, aes(x = est, fill = sig_fac)) +
  geom_histogram(breaks = seq(-1*max_slope, max_slope, by = 0.1),
                 position = "stack",
                 colour = "white", show.legend = FALSE) +
  scale_fill_manual(values = c(alpha("steelblue3", 0.8), alpha("steelblue3", 0.3))) +
  theme_bw() +
  geom_vline(xintercept = 0) +
  ylab("Number of species") +
  xlab("Linear model slope") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(path(figures_dir, "figure-03a_bam-linear-slope_histogram.png"),
       width = 14, height = 10, units = "cm")
ggsave(path(figures_dir, "figure-03a_bam-linear-slope_histogram.tif"),
       width = 14, height = 10, units = "cm")


## summarise significance and direction
spec_coef <- spec_coef |> 
  mutate(direction = ifelse(est>0, "pos", "neg"))
spec_coef |>
  dplyr::select(sig_fac, direction) |>
  table()
#        direction
# sig_fac neg pos
#      s  359  49
#      ns  50  37

sum(spec_coef$direction == "neg") / nrow(spec_coef)
# [1] 0.8279352
sum(spec_coef$direction == "neg" & spec_coef$sig_fac == "s") / nrow(spec_coef)
# [1] 0.7267206




#########################################################
## biome summaries

biome_coef |> arrange(est)

#   breeding_biome     int_est   est min_log_abd max_log_abd
#   <chr>                <dbl> <dbl>       <dbl>       <dbl>
# 1 Aridlands           2.43   -2.72           0           1
# 2 Grassland           0.0436 -2.47           0           1
# 3 Wetlands & Coasts  -0.0557 -2.23           0           1
# 4 Arctic Tundra      -1.39   -1.93           0           1
# 5 Habitat Generalist  0.579  -1.66           0           1
# 6 Forest              0.784  -1.64           0           1


spec_coef |> 
  mutate(direction = ifelse(est < 0, "neg", "pos")) |>
  mutate(sig_neg = ifelse(direction == "neg" & sig == "s", 1, 0)) |>
  mutate(neg_bin = ifelse(est < 0, 1, 0)) |>
  dplyr::select(breeding_biome, sig_neg, neg_bin) |>
  group_by(breeding_biome) |>
  summarise_all(mean) |>
  arrange(neg_bin)


# proportion of species in each biome with negative slopes 
# and with significantly negative slopes

#   breeding_biome     sig_neg neg_bin
#   <chr>                <dbl>   <dbl>
# 1 Arctic Tundra        0.7     0.733
# 2 Habitat Generalist   0.780   0.805
# 3 Forest               0.700   0.808
# 4 Wetlands & Coasts    0.680   0.836
# 5 Aridlands            0.789   0.868
# 6 Grassland            0.913   0.957


#########################################################
## 
# figure 10-year cumulative population change with abundance, 
# plotted by biome. 


nu <- function(x){as.numeric(as.character(x))}
plot_1spec_line_multiyr <- function(data, col = alpha("black", 0.1), no_years = 10, plot = TRUE, ...){
  x <- seq(0, 1, by = 0.004)
  x2 <- seq(data$min_log_abd[1], data$max_log_abd[1], length.out = length(x))
  tr_1yr <- x2*data$est + data$int_est
  tr_multiyr <- ppy_to_multiyr(tr_1yr, noyrs = no_years)
  if(plot) lines(x, tr_multiyr, col = col, ...)
  if(!plot) return(data.frame(x, tr_multiyr))
}

plot_multispec_multiyr <- function(data, col = alpha("black", 0.1), no_years = 10){
  for(i in 1:nrow(data)){
    plot_1spec_line_multiyr(data[i,], col = col, no_years = no_years)
  }
}


no_yrs <- 10
ylim <- 60
ytick <- 50


# change breeding biome names to match paper
biomes <- names(table(spec_coef$breeding_biome))
biome_names_tidy <- c("Arctic tundra", "Aridland", "Forest",
                      "Grassland", "Habitat generalists", "Wetland & Coast")

png(path(figures_dir, "figure-03b_bam-linear-slopes.png"),
    width = 14.4, height = 10, units = "cm", pointsize = 9, res = 600)

par(mfrow=c(2, 3), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,6,2,0))


for(i in 1:length(biomes)) {
  ylabel <- paste("Population change over", no_yrs, "years (percent)")
  plot(0, 0, col = "white", ylim = ylim*c(-1, 1), xlim = c(0, 1),
       main = "", xaxt = "n", xlab = "", yaxt = "n", ylab = "",
       xpd = NA, cex.lab = 1.5, las = 1)
  if(i %in% c(1,4)) axis(side = 2, at = ytick * c(-1, 0, 1), cex.axis = 1.5, las = 1)
  if(i == 4) text(x = -0.3, y = 1.1*ylim, label = ylabel, cex = 1.5, xpd = NA, srt = 90)
  if(i > 3) axis(side = 1, at = c(0, 1), labels = c("0", "max"), cex.axis = 1.5, las = 1)
  abline(h = 0, col = alpha("firebrick", 0.42), lwd = 2)
  
  # lines for each species in the biome
  plot_multispec_multiyr(filter(spec_coef, breeding_biome == biomes[i]), no_years = no_yrs)
  
  if(i == 5) text(x = 0.5, y = -1.66*ylim,
                  labels = "Relative abundance within species", xpd = NA, cex = 1.5)
  text(x = -0.03, y = 11/12*ylim, pos = 4, labels = biome_names_tidy[i], cex = 1.5)
  
  # add average biome line
  biome_avg <- biome_coef |> filter(breeding_biome == biomes[i])
  plot_1spec_line_multiyr(biome_avg, col = "black", lwd = 2, no_years = no_yrs)
}

dev.off()


diff10 <- vector(length = length(biomes))

for(i in 1:6){
  biome_avg <- biome_coef |> filter(breeding_biome == biomes[i])
  pred10 <- plot_1spec_line_multiyr(biome_avg, col = "black", lwd = 2, no_years = no_yrs, plot = FALSE)
  diff10[i] <- pred10$tr_multiyr[nrow(pred10)] - pred10$tr_multiyr[1]
}

cbind(biomes, diff10)


#########################################################
## individual species plots

focal_species <- c("wilsap", "yebcuc")

species_names <- auk::get_ebird_taxonomy(2022) |> 
  select(species_code, common_name) |>
  filter(species_code %in% focal_species)

for(i in 1:nrow(species_names)){
  
  # collect species data for plots
  this_species_code <- species_names$species_code[i]
  this_species_name <- species_names$common_name[i]
  this_species_coef <- spec_coef |> filter(species_code == this_species_code)
  this_species_trends <- trends |> filter(species_code == this_species_code)
  
  # set up plot (to hopefully make pretty for many species)
  xlimits <- range(this_species_trends$log_abd)
  la_vals <- seq(-5, 5, by = 1)
  
  #  if(xlimits[2] - xlimits[1])
  #  la_vals <- seq(-6, 6, by = 2)
  
  ylimits <- range(this_species_trends$abd_ppy_median)
  if(ylimits[2]< -10) ylimits[2] <- -10
  if(ylimits[1]> 10) ylimits[1] <- 10
  ylimits <- c(-10, 10)
  
  
  formatted_slope <- paste0(ifelse(nu(this_species_coef$est)<0, "n", "p"), "_",
                            format(round(abs(nu(this_species_coef$est)), digits = 3), nsmall = 3))
  
  plot_loc <- path(figures_dir, paste0("species_scatter_", mod_tag, "_",
                                       this_species_code, ".png"))
  png(plot_loc, width = 8, height = 8, units = "cm", pointsize = 9, res = 600)
  par(mar = c(5, 5, 3, 3))
  
  # create plot
  plot(0, 0, col = "white", 
       xlim = xlimits, 
       xaxt = "n", xlab = "Relative abundance",
       ylim = ylimits,
       yaxt = "n", ylab = "Population trend", 
       main = this_species_name, 
       mgp = c(1.5, 0.7, 0))
  axis(side = 1, at = la_vals, labels = 10^la_vals, tcl = -0.3, mgp = c(1.5, 0.5, 0))
  axis(side = 2, at = c(-10, 0, 10), las = 1, tcl = -0.3, mgp = c(1.5, 0.7, 0))
  
  abline(h = 0, col = alpha("firebrick", 0.2), lwd = 2)
  
  points(nu(this_species_trends$log_abd), nu(this_species_trends$abd_ppy_median), 
         pch = 16, col = alpha("steelblue", 0.3),
         cex = sqrt(nu(this_species_trends$abd_ppy_sd / mean(this_species_trends$abd_ppy_sd) )))
  
  abline(nu(this_species_coef$int_est), nu(this_species_coef$est), col = "steelblue", lwd= 1.5)
  
  dev.off()
}
