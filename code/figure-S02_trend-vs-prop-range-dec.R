library(tidyverse)
library(lme4)
library(arrow)
library(fs)

data_dir <- "data"
outputs_dir <- "outputs"
figures_dir <- path("figures", "figure-S02")
dir_create(figures_dir)


# load data ----

# trends estimates
trends <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet() |>
  mutate(log_abd = log10(abd),
         breeding_biome = factor(breeding_biome),
         species_code = factor(species_code)) |>
  dplyr::select(species_code, season, abd, abd_ppy_median)

# summarise by species: proportion of declining cells and rangewide trend
by_species <- trends |>
  mutate(abd_mult_ppy = abd * abd_ppy_median) |>
  group_by(species_code, season) |>
  summarise(tot_abd = sum(abd), 
            tot_wt_ppy = sum(abd_mult_ppy),
            prop_dec = mean(abd_ppy_median < 0),
            .groups = "drop") |>
  mutate(wt_trend = tot_wt_ppy / tot_abd) |>
  mutate(wt_trend_cat = cut(wt_trend, breaks = seq(-10, 4, by = 1)))


# boxplot ----

plot_loc <- path(figures_dir, paste0("figure-S02_boxplot.tif"))
tiff(plot_loc, width = 8, height = 8, units = "cm", pointsize = 9, res = 600)
par(mar = c(5, 5, 3, 3))

boxplot(prop_dec ~ wt_trend_cat, data = by_species, 
        range = 0, lty = 1, xaxt = "n", yaxt = "n", 
        xlab = "Range wide trend", ylab = "Proportion of range with declines")
axis(side = 2, at = c(0, 0.5, 1), ylab = "Proportion of range with declines")
axis(side = 1, at = seq(0.5, 14.5, by = 1), labels = c(-10:4))
abline(v = 10.5, col = "red")

dev.off()
