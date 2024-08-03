
library(tidyverse)
library(arrow)
library(fs)

base_dir <- "/Users/Alison/Documents/REPOS/2024_trends-abundance-paper/"
data_dir <- path(base_dir, "data")
outputs_dir <- path(base_dir, "outputs")
figures_dir <- path(base_dir, "figures/figure_S1")
dir.create(figures_dir)


#########################################################
## load data

# bcr lookup
bcr_lookup <- path(data_dir, "bcr-srd-lookup.csv") |>
  read_csv() |>
  filter(coverage_fraction > 0.5)

# trends estimates
trends <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet() |>
  inner_join(bcr_lookup, by = join_by("srd_id")) |>
  dplyr::select(species_code, abd_ppy_median, region_code)


#########################################################
## summarise by species and bcr

tr_sp_bcr <- trends |>
  group_by(species_code, region_code) |>
  summarise(tr_min = min(abd_ppy_median), tr_max = max(abd_ppy_median), no_cells = n()) |>
  mutate(tr_range = tr_max - tr_min) |>
  mutate(at_least_20 = ifelse(no_cells >= 20, 1, 0))

med_all <- median(tr_sp_bcr$tr_range)
med_20 <- median(tr_sp_bcr$tr_range[tr_sp_bcr$at_least_20 == 1])


#########################################################
## plot histogram of range of species trends within bcr 

h1 <- hist(tr_sp_bcr$tr_range, breaks = 30, plot = FALSE)

plot_loc <- path(figures_dir, "figure_S1_histogram_species_bcr_ranges.tif")
tiff(plot_loc, width = 12, height = 8, units = "cm", pointsize = 9, res = 600)
par(mar = c(5, 5, 1, 1))

  hist(tr_sp_bcr$tr_range, breaks = h1$breaks, col = "grey70", 
    xlim = c(0, 40), xlab = "Range of single-species trends within each BCR (percent-per-year)", xaxt = "n",
    yaxt = "n", main = "")
  axis(side = 1, at = c(0, 10, 20, 30, 40))
  axis(side = 2, at = c(0, 1000, 2000))
  hist(tr_sp_bcr$tr_range[tr_sp_bcr$at_least_20==1], breaks = h1$breaks, col = "grey30", add = TRUE)
  segments(x0 = med_all, x1 = med_all, y0 = 0, y1 = 3000, col = "white", lty = 2)

  legend(x = 20, y = 2000, fill = c("grey70", "grey30"), 
    legend = c("BCR-species < 20 cells", "BCR-species >= 20 cells"), cex = 0.8)

dev.off()







