library(tidyverse)
library(arrow)
library(fs)

data_dir <-  "data"

# trends estimates
trends <- path(data_dir, "ebird-trends_2007-2021.parquet") |>
  read_parquet() |>
  dplyr::select(species_code, abd, abd_ppy_median)

# summarise by species and calculate weighted range-wide trend
tr_sp <- trends |>
  group_by(species_code) |>
  summarise(tr_min = min(abd_ppy_median), tr_max = max(abd_ppy_median), no_cells = n(),
    tr_med = median(abd_ppy_median), sum_abd = sum(abd), sum_tr_abd = sum(abd*abd_ppy_median),
    prop_range_dec = mean(abd_ppy_median < 0),
    .groups = "drop") |>
  mutate(tr_range = tr_max - tr_min,
          tr_wt_avg = sum_tr_abd / sum_abd) |>
  dplyr::select(species_code, no_cells, tr_med, tr_wt_avg, tr_range, prop_range_dec) |>
  arrange(tr_wt_avg)

mean(tr_sp$tr_wt_avg<0) # 0.7515152
mean(tr_sp$tr_med<0)    # 0.6747475
mean(tr_sp$prop_range_dec) # 0.6116947

tr_sp |> filter(species_code %in% c("rufhum", "balori"))
