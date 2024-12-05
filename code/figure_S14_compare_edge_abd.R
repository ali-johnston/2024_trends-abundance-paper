#########################################################
## fit the linear mixed models for abundance and edge
## create figure S13

library(tidyverse)
library(fs)
library(lme4)
library(arrow)
library(performance)



base_dir <- "/Users/Alison/Documents/REPOS/2024_trends-abundance-paper/"
data_dir <- path(base_dir, "data")
outputs_dir <- path(base_dir, "outputs")
figures_dir <- path(base_dir, "figures/figure_S14")
dir.create(figures_dir)



#########################################################
## load data

# # species lookup
species <- path(data_dir, "master_species_list_495.csv") |>
  read_csv(na = "") |>
  select(species_code, breeding_biome)


# trends estimates
trends <- path(data_dir, "ebird-trends_2021_weights.parquet") |>
  read_parquet() |>
  left_join(species) |>
  mutate(log10_abd = log10(abd),
         log10_distance_to_edge_km = log10(distance_to_edge_km),
         breeding_biome = factor(breeding_biome),
         species_code = factor(species_code),
         weight = abd_ppy_sd / mean(abd_ppy_sd))


cor(trends$log10_abd, trends$log10_distance_to_edge_km)
# [1] 0.5173799

dist_cor <- trends |> 
  dplyr::select(species_code, log10_abd, log10_distance_to_edge_km) |>
  group_by(species_code) |>
  summarise(adc = cor(log10_abd, log10_distance_to_edge_km)) |>
  arrange(adc) |>
  mutate(adc2 = adc^2)

median(dist_cor$adc2)
# [1] 0.2869478


#########################################################
## calculate min and max abundance

# minimum and maximum abundance for each species and each biome
abd_range <- trends |>
  dplyr::select(species_code, breeding_biome, abd) |>
  group_by(species_code, breeding_biome) |>
  summarise(max_abd = max(abd), min_abd = min(abd), .groups = "drop") |>
  mutate(min_log_abd = log10(min_abd), max_log_abd = log10(max_abd))


#########################################################
## read in abundance results and calculate effect sizes

mod_tag <- "log10_abd_40"

results_loc1 <- path(outputs_dir, paste0("species_coefs_bam_linear_", mod_tag, ".csv"))
results_loc2 <- path(outputs_dir, paste0("species_coefs_lm_linear_", mod_tag, ".csv"))

spec_coef_abd <- read_csv(results_loc1) |>
              bind_rows(read_csv(results_loc2)) |>
              mutate(sig = ifelse(p_val < 0.05, "s", "ns")) |>
              mutate(sig_fac = factor(sig, levels = c("s", "ns"), ordered = TRUE)) |>
              rename(species_code = species_code...5) |>
              left_join(species, by = "species_code") |>
              dplyr::select(-species_code...10) |>
              merge(dplyr::select(abd_range, -breeding_biome), by = "species_code") |>
              mutate(effect_size = est*(max_log_abd - min_log_abd)) |>
              dplyr::select(species_code, effect_size) |>
              rename(pred_diff_abd = effect_size)


#########################################################
## calculate min and max distance to range edge

dist_range <- trends |>
  dplyr::select(species_code, breeding_biome, distance_to_edge_km) |>
  group_by(species_code, breeding_biome) |>
  summarise(max_dist = max(distance_to_edge_km), min_dist = min(distance_to_edge_km), .groups = "drop") |>
  mutate(min_log_dist = log10(min_dist), max_log_dist = log10(max_dist))


#########################################################
## read in abundance results and calculate effect sizes

mod_tag <- "log10_distance_to_edge_km_40"

results_loc1 <- path(outputs_dir, paste0("species_coefs_bam_linear_", mod_tag, ".csv"))
results_loc2 <- path(outputs_dir, paste0("species_coefs_lm_linear_", mod_tag, ".csv"))

spec_coef_edge <- read_csv(results_loc1) |>
              bind_rows(read_csv(results_loc2)) |>
              mutate(sig = ifelse(p_val < 0.05, "s", "ns")) |>
              mutate(sig_fac = factor(sig, levels = c("s", "ns"), ordered = TRUE)) |>
              rename(species_code = species_code...5) |>
              left_join(species, by = "species_code") |>
              dplyr::select(-species_code...10) |>
              merge(dplyr::select(dist_range, -breeding_biome), by = "species_code") |>
              mutate(effect_size = est*(max_log_dist - min_log_dist)) |>
              dplyr::select(species_code, effect_size) |>
              rename(pred_diff_edge = effect_size)



#########################################################
## compare diffs in min to max predictions for two variables

pred_compare <- spec_coef_abd |>
      inner_join(spec_coef_edge)


rng <- range(pred_compare[,c("pred_diff_abd", "pred_diff_edge")])

br <- seq(rng[1], rng[2], length.out = 30)
h1 <- hist(pred_compare$pred_diff_edge, breaks = br, plot = FALSE)
h2 <- hist(pred_compare$pred_diff_abd, breaks = br, plot = FALSE)
mc <- max(h1$counts, h2$counts)



plot_loc <- path(figures_dir, "figure_S14_effect_size_edge_abd.tif")
plot_loc <- path(figures_dir, "figure_S14_effect_size_edge_abd.png")
png(plot_loc, width = 12, height = 12, units = "cm", pointsize = 9, res = 600)

  par(mfrow = c(2,1), mar = c(1, 5, 0.5, 0.5), oma = c(4, 0, 0, 0))

  hist(pred_compare$pred_diff_edge, breaks = br, main = "", xlab = "", xaxt = "n", ylim = c(0, mc))
  text(x = rng[1], 
    y = mc*0.9, 
    labels = "a) Distance to edge (km)", font = 2, pos = 4)
  axis(side = 1, at = c(-20, -10, 0, 10), labels = rep("", 4))

  hist(pred_compare$pred_diff_abd, breaks = br, main = "", xlab = "Effect size of variable on trend", xaxt = "n", xpd = NA, ylim = c(0, mc))
  text(x = rng[1], 
    y = mc*0.9, 
    labels = "b) Relative abundance", font = 2, pos = 4)
  axis(side = 1, at = c(-20, -10, 0, 10))

dev.off()



#########################################################
## testing differences between distributions. 

t.test(pred_compare$pred_diff_edge, mu = 0, alternative = "two.sided")
#   One Sample t-test

# data:  pred_compare$pred_diff_edge
# t = -2.7791, df = 494, p-value = 0.005659
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  -0.25569871 -0.04389299
# sample estimates:
#  mean of x 
# -0.1497959 

t.test(pred_compare$pred_diff_abd, mu = 0, alternative = "two.sided")
#   One Sample t-test

# data:  pred_compare$pred_diff_abd
# t = -17.684, df = 494, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  -2.232578 -1.786095
# sample estimates:
# mean of x 
# -2.009337 


t.test(pred_compare$pred_diff_edge, pred_compare$pred_diff_abd)
#   Welch Two Sample t-test

# data:  pred_compare$pred_diff_edge and pred_compare$pred_diff_abd
# t = 14.787, df = 705.62, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  1.612635 2.106446
# sample estimates:
#  mean of x  mean of y 
# -0.1497959 -2.0093365 

library(car)

df <- pred_compare |>
      dplyr::select(pred_diff_abd, pred_diff_edge) |>
      pivot_longer(cols = 1:2, names_to = "type", values_to = "pred_diff") |>
      mutate(type = ifelse(type=="pred_diff_abd", "abd", "edge"))

leveneTest(pred_diff ~ type, data = df)
# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
# group   1  138.21 < 2.2e-16 ***
#       988                      
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Warning message:
# In leveneTest.default(y = y, group = group, ...) : group coerced to factor.

sd(pred_compare$pred_diff_edge)
# [1] 1.199214
sd(pred_compare$pred_diff_abd)
# [1] 2.527926


