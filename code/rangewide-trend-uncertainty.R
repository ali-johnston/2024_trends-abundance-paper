library(arrow)
library(dplyr)
library(fs)
library(purrr)
library(readr)

data_dir <- "data"

# read in trends estimates with breeding biomes and srd information
rw_folds <- path(data_dir, "ebird-trends_range-wide_folds_2021.parquet") |>
  read_parquet()


# calculate uncertainty based on folds 
q20 <- function(x) {
  quantile(x$rw_ppy, probs = c(0.1, 0.5, 0.9))
}
rw_uncert <- rw_folds |>
  dplyr::select(species_code, season, rw_ppy = abd_ppy) |>
  group_by(species_code, season) |>
  group_split() |>
  map_dfr(q20) |>
  rename(q10 = '10%', q50 = '50%', q90 = '90%') |>
  mutate(dir = ifelse(q10>0, "sig_pos",
                      ifelse(q90<0, "sig_neg",
                             ifelse(q50<0, "neg", "pos")))) |>
  pull(dir) |>
  table() / 495