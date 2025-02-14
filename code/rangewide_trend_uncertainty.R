library(dplyr)
library(fs)
library(tidyverse)

base_dir <- "/Users/Alison/Documents/REPOS/2024_trends-abundance-paper/"
data_dir <- path(base_dir, "data")
outputs_dir <- path(base_dir, "outputs")
figures_dir <- path(base_dir, "figures")



#########################################################
## read in data

# # species lookup
species <- read_csv(path(data_dir, "master_species_list_495.csv"), na = "") |>
  select(species_code, season, breeding_biome) |>
  mutate(run_name = paste(species_code, season, sep = "_"))


# read in trends estimates with breeding biomes and srd information
rw_folds <- path(data_dir, "rangewide-trends_folds.csv") |>
  read_csv() |>
  right_join(species, by = "run_name")


#########################################################
## calculate uncertainty based on folds 

q20 <- function(x) {quantile(x$rw_ppy, probs = c(0.1, 0.5, 0.9))}

rw_uncert <- rw_folds |>
        dplyr::select(run_name, rw_ppy) |>
        group_by(run_name) |>
        group_split() |>
        map_dfr(q20) |>
        rename(q10 = '10%', q50 = '50%', q90 = '90%') |>
        mutate(dir = ifelse(q10>0, "sig_pos", ifelse(q90<0, "sig_neg", ifelse(q50<0, "neg", "pos")))) |>
        pull(dir) |>
        table() / 495

