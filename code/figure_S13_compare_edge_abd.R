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
figures_dir <- path(base_dir, "figures/figure_S13")
dir.create(figures_dir)



#########################################################
## load data

# # species lookup
species <- path(data_dir, "master_species_list_495.csv") |>
  read_csv(na = "") |>
  select(species_code, breeding_biome)


# trends estimates
trends <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet() |>
  mutate(log_abd = log10(abd),
         log_distance_to_edge_km = log10(distance_to_edge_km),
         breeding_biome = factor(breeding_biome),
         species_code = factor(species_code),
         weight = 1 / abd_ppy_var)


cor(trends$log_abd, trends$log_distance_to_edge_km)

dist_cor <- trends |> 
  dplyr::select(species_code, log_abd, log_distance_to_edge_km) |>
  group_by(species_code) |>
  summarise(adc = cor(log_abd, log_distance_to_edge_km)) |>
  arrange(adc) |>
  mutate(adc2 = adc^2)




#########################################################
## Mixed model with abundance and 1/var as weight: 

# mixed model
m_abd <- lmer(abd_ppy_median ~ log_abd + 
          (1 + log_abd|breeding_biome/species_code),
           weights = weight,
           data = trends)

m_edge <- lmer(abd_ppy_median ~ log_distance_to_edge_km + 
          (1 + log_distance_to_edge_km|breeding_biome/species_code),
           weights = weight,
           data = trends)

# the slope of log_distance_to_edge_km on breeding_biome
# had a low variance near to zero. 
# subsequent model performance checks did not work. 

# so changing the model structure to remove this: 

m_edge2 <- lmer(abd_ppy_median ~ log_distance_to_edge_km + 
          (1 + log_distance_to_edge_km|species_code),
           weights = weight,
           data = trends)

# remove the breeding biome from the abundance model to 
# allow a like-for-like comparison

m_abd2 <- lmer(abd_ppy_median ~ log_abd + 
          (1 + log_abd|species_code),
           weights = weight,
           data = trends)


AIC(m_edge2) - AIC(m_abd2)
# [1] 67494.1

#########################################################
## minimum and maximum of the distance to edge variable for each species and each biome
## then predictions for the max and min

mod <- m_edge2

edge_df0 <- trends |>
  dplyr::select(species_code, breeding_biome, log_distance_to_edge_km, log_abd) |>
  group_by(species_code, breeding_biome) |>
  summarise(max_edge = max(log_distance_to_edge_km), min_edge = min(log_distance_to_edge_km), 
    q0.25_edge = quantile(log_distance_to_edge_km, probs = 0.25), q0.75_edge = quantile(log_distance_to_edge_km, probs = 0.75), 
    med_abd = median(log_abd), .groups = "drop") |>
  pivot_longer(cols = ends_with("_edge"),
            names_to = "minmax",
            values_to = "log_distance_to_edge_km") |>
  rename(log_abd = med_abd)

edge_df0$pred <- predict(mod, newdata = edge_df0)

edge_pred <- edge_df0 |>
  dplyr::select(species_code, breeding_biome, minmax, pred) |>
  pivot_wider(names_from = minmax, 
            values_from = pred) |>
  mutate(pred_diff_edge = max_edge - min_edge) |>
  mutate(pred_diff_edge_iqr = q0.75_edge - q0.25_edge) |>
  dplyr::select(species_code, breeding_biome, pred_diff_edge, pred_diff_edge_iqr)



#########################################################
## minimum and maximum of the abundance variable for each species and each biome
## then predictions for the max and min

mod <- m_abd2

abd_df0 <- trends |>
  dplyr::select(species_code, breeding_biome, log_distance_to_edge_km, log_abd) |>
  group_by(species_code, breeding_biome) |>
  summarise(max_abd = max(log_abd), min_abd = min(log_abd), 
    q0.25_abd = quantile(log_abd, probs = 0.25), q0.75_abd = quantile(log_abd, probs = 0.75),
    med_edge = median(log_distance_to_edge_km), .groups = "drop") |>
  pivot_longer(cols = ends_with("_abd"),
            names_to = "minmax",
            values_to = "log_abd") |>
  rename(log_distance_to_edge_km = med_edge)

abd_df0$pred <- predict(mod, newdata = abd_df0)

abd_pred <- abd_df0 |>
  dplyr::select(species_code, breeding_biome, minmax, pred) |>
  pivot_wider(names_from = minmax, 
            values_from = pred) |>
  mutate(pred_diff_abd = max_abd - min_abd) |>
  mutate(pred_diff_abd_iqr = q0.75_abd - q0.25_abd) |>
  dplyr::select(species_code, breeding_biome, pred_diff_abd, pred_diff_abd_iqr)


#########################################################
## compare diffs in min to max predictions for two variables

pred_compare <- abd_pred |>
      inner_join(edge_pred)


rng <- range(pred_compare[,c("pred_diff_abd", "pred_diff_edge")])

br <- seq(rng[1], rng[2], length.out = 30)
h1 <- hist(pred_compare$pred_diff_edge, breaks = br, plot = FALSE)
h2 <- hist(pred_compare$pred_diff_abd, breaks = br, plot = FALSE)
mc <- max(h1$counts, h2$counts)



plot_loc <- path(figures_dir, "figure_S13_effect_size_edge_abd.tif"))
tiff(plot_loc, width = 12, height = 12, units = "cm", pointsize = 9, res = 600)

  par(mfrow = c(2,1), mar = c(1, 5, 0.5, 0.5), oma = c(4, 0, 0, 0))

  hist(pred_compare$pred_diff_edge, breaks = br, main = "", xlab = "", xaxt = "n", ylim = c(0, mc))
  text(x = rng[1], 
    y = mc*0.9, 
    labels = "a) Distance to edge (km)", font = 2, pos = 4)
  axis(side = 1, at = c(-20, -10, 0, 10), labels = rep("", 4))

  hist(pred_compare$pred_diff_abd, breaks = br, main = "", xlab = "Difference in species trend explained by variable", xaxt = "n", xpd = NA, ylim = c(0, mc))
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
# t = -1.5937, df = 494, p-value = 0.1116
# alternative hypothesis: true mean is not equal to 0

t.test(pred_compare$pred_diff_abd, mu = 0, alternative = "two.sided")
#   One Sample t-test

# data:  pred_compare$pred_diff_abd
# t = -13.393, df = 494, p-value < 2.2e-16

t.test(pred_compare$pred_diff_edge, pred_compare$pred_diff_abd)
#   Welch Two Sample t-test

# data:  pred_compare$pred_diff_edge and pred_compare$pred_diff_abd
# t = 11.071, df = 753.49, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0

library(car)

df <- pred_compare |>
      dplyr::select(pred_diff_abd, pred_diff_edge) |>
      pivot_longer(cols = 1:2, names_to = "type", values_to = "pred_diff") |>
      mutate(type = ifelse(type=="pred_diff_abd", "abd", "edge"))

leveneTest(pred_diff ~ type, data = df)
# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
# group   1  95.385 < 2.2e-16 ***
#       988                      
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

sd(pred_compare$pred_diff_edge)
# [1] 1.573363
sd(pred_compare$pred_diff_abd)
# [1] 2.953419


