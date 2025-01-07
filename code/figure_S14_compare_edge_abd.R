#########################################################
## fit the linear mixed models for abundance and edge
## create figure S13

library(tidyverse)
library(fs)
library(lme4)
library(arrow)
library(performance)
library(car)

data_dir <- "data"
outputs_dir <- "outputs"
figures_dir <- path("figures", "figure_S14")
dir_create(figures_dir)


#########################################################
## load data

# # species lookup
species <- path(data_dir, "master_species_list_495.csv") |>
  read_csv(na = "") |>
  select(species_code, breeding_biome)


# trends estimates
trends <- path(data_dir, "ebird-trends_2021_weights.parquet") |>
  read_parquet() |>
  left_join(species, by = join_by(species_code)) |>
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

mod_tag <- "log10_abd_40_run2"

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

mod_tag <- "log10_distance_to_edge_km_40_run2"

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




#########################################################
## compare AIC


# -------------------------------------------------------
# read in abundance model AIC

mod_tag <- "log10_abd_40_run2"

results_loc1_abd <- path(outputs_dir, paste0("species_coefs_bam_linear_", mod_tag, ".csv"))
results_loc2_abd <- path(outputs_dir, paste0("species_coefs_lm_linear_", mod_tag, ".csv"))

aic_abd_df <- read_csv(results_loc1_abd) |>
              bind_rows(read_csv(results_loc2_abd)) |>
              dplyr::select(species_code...5, aic) |>
              rename(species_code = species_code...5) |>
              rename(aic_abd = aic)

# -------------------------------------------------------
# read in edge model AIC

mod_tag <- "log10_distance_to_edge_km_40_run2"

results_loc1_edge <- path(outputs_dir, paste0("species_coefs_bam_linear_", mod_tag, ".csv"))
results_loc2_edge <- path(outputs_dir, paste0("species_coefs_lm_linear_", mod_tag, ".csv"))

aic_edge_df <- read_csv(results_loc1_edge) |>
              bind_rows(read_csv(results_loc2_edge)) |>
              dplyr::select(species_code...5, aic) |>
              rename(species_code = species_code...5) |>
              rename(aic_edge = aic)


# -------------------------------------------------------
# compare AIC by species

aic_comp <- aic_abd_df |>
              left_join(aic_edge_df, by = "species_code") |>
              mutate(delta_aic = aic_edge - aic_abd) |>
              mutate(abd_imp = ifelse(delta_aic > 0, 1, 0)) |>
              mutate(delta_aic_within2 = ifelse(delta_aic < 2 & delta_aic > -2, 1, 0)) |>
              mutate(abd_imp_2 = ifelse(delta_aic > 2, 1, 0)) |>
              mutate(edge_imp_2 = ifelse(delta_aic < -2, 1, 0)) |>
              filter(!is.na(aic_abd), !is.na(aic_edge))

nrow(aic_comp)

median(aic_comp$delta_aic, na.rm = TRUE)
# [1] 28.15127
sd(aic_comp$delta_aic, na.rm = TRUE)
# [1] 254.7087

mean(aic_comp$delta_aic_within2, na.rm = TRUE) |> round(2)
# [1] 0.09

mean(aic_comp$abd_imp_2, na.rm = TRUE) |> round(2)
# [1] 0.71

mean(aic_comp$edge_imp_2, na.rm = TRUE) |> round(2)
# [1] 0.20

hist(aic_comp$delta_aic, breaks = 100)
hist(aic_comp$delta_aic, breaks = seq(-2000, 3000, by = 2), xlim = c(-50, 50))


#########################################################
## AIC by effect size

aic_effect_size <- aic_comp |>
            dplyr::select(species_code, delta_aic) |>
            left_join(spec_coef_abd, by = "species_code") |>
            mutate(neg = ifelse(pred_diff_abd < 0, 1, 0)) |>
            mutate(delta_cat = ifelse(delta_aic < -2, "A", ifelse(delta_aic < 2, "B", "C"))) |>
            mutate(delta_cat2 = ifelse(delta_aic < -4, "A", ifelse(delta_aic < -2, "B", ifelse(delta_aic < 2, "C", ifelse(delta_aic < 4, "D", "E")))))

aic_effect_size |>
          dplyr::select(delta_cat, neg) |>
          group_by(delta_cat) |> 
          summarise(prop_neg = mean(neg), n_spec = n())

# # A tibble: 3 x 3
#   delta_cat prop_neg n_spec
#   <chr>        <dbl>  <int>
# 1 A            0.526     97  0.53
# 2 B            0.766     47  0.77
# 3 C            0.917    351  0.92


#########################################################
## read in abundance results and calculate effect sizes

mod_tag <- "log10_abd_40_run2"

results_loc1 <- path(outputs_dir, paste0("species_coefs_bam_linear_", mod_tag, ".csv"))
results_loc2 <- path(outputs_dir, paste0("species_coefs_lm_linear_", mod_tag, ".csv"))

aic_sig_cross <- read_csv(results_loc1) |>
              bind_rows(read_csv(results_loc2)) |>
              mutate(sig = ifelse(p_val < 0.05, "s", "ns")) |>
              mutate(sig_fac = factor(sig, levels = c("s", "ns"), ordered = TRUE)) |>
              rename(species_code = species_code...5) |>
              left_join(species, by = "species_code") |>
              dplyr::select(-species_code...10) |>
              merge(dplyr::select(abd_range, -breeding_biome), by = "species_code") |>
              mutate(effect_size = est*(max_log_abd - min_log_abd)) |>
              dplyr::select(species_code, effect_size, sig, sig_fac) |>
              rename(pred_diff_abd = effect_size) |>
              left_join(dplyr::select(aic_comp, species_code, delta_aic, abd_imp, delta_aic_within2)) |>
              mutate(delta_cat = ifelse(delta_aic < -2, "edge", ifelse(delta_aic < 2, "both", "abd"))) |>
              mutate(dir = ifelse(pred_diff_abd < 0, ifelse(sig == "s", "neg_sig", "neg"), ifelse(sig == "s", "pos_sig", "pos"))) |>
              mutate(dir = factor(dir, levels = c("neg_sig", "neg", "pos", "pos_sig"), ordered = TRUE)) |>
              dplyr::select(delta_cat, dir) |>
              table()

aic_sig_cross

#          dir
# delta_cat neg_sig neg pos pos_sig
#      abd      171  78  62      40
#      both       7  15  19       6
#      edge      17   2   3      75

chisq.test(aic_sig_cross)
#   Pearson's Chi-squared test

# data:  aic_sig_cross
# X-squared = 211.51, df = 6, p-value < 2.2e-16


rowsums <- matrix(rep(apply(aic_sig_cross, 1, sum) /495, 4), 
  nrow = nrow(aic_sig_cross), ncol = ncol(aic_sig_cross))

colsums <- matrix(rep(apply(aic_sig_cross, 2, sum) /495, 3), 
  nrow = nrow(aic_sig_cross), ncol = ncol(aic_sig_cross), byrow = TRUE)

expected <- round(rowsums * colsums * 495, 3)

o_minus_e <- round(aic_sig_cross - expected, 0)

# delta_cat neg_sig neg pos pos_sig
#      abd       33  11   2     -46
#      both     -12   6  11      -5
#      edge     -21 -17 -13      51



