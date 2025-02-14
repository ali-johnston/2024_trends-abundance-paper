library(arrow)
library(dplyr)
library(fs)
library(ggplot2)
library(purrr)
library(tidyr)

data_dir <- "data"
figures_dir <- path("figures", "figure-S09-S10")
dir_create(figures_dir)


# load data ----

# directional error,focal species, trend bins
abd <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet() |>
  distinct(species_code, season, srd_id, abd)
trends <- path(data_dir, "ebird-trends_simulations_focal-species_2021.parquet") |>
  read_parquet() |>
  inner_join(abd, by = join_by(species_code, season, srd_id)) |> 
  # assign trend direction
  mutate(direction_sim = sign(simulated),
         direction_est = sign(estimated),
         correct_direction = nonzero & direction_sim == direction_est,
         incorrect_direction = nonzero & direction_sim != direction_est)


# model performance function ----

model_performance <- function(data) {
  s <- data.frame(n_nonzero = sum(data$nonzero),
                  directional_power = mean(data$correct_direction))
  data_nz <- data[data$nonzero, ]
  if (nrow(data_nz) > 0) {
    s$directional_error <- mean(data_nz$incorrect_direction)
  } else {
    s$directional_error <- NA_real_
  }
  return(s)
}


# calculate dinned model performance
species <- unique(trends$species_code)
tax <- auk::ebird_taxonomy
ppy_binned <- NULL
abd_binned <- NULL
for (s in species) {
  common_name <- tax$common_name[tax$species_code == s]
  t <- filter(trends, species_code == s)
  
  # ppy magnitude bins, 0.5 width
  t$ppy_bin_midpoint <- 0.5 * (abs(t$estimated) %/% 0.5) + 0.5 / 2
  ppy_binned <- t |>
    nest(data = !ppy_bin_midpoint) |>
    mutate(performance = map(data, model_performance)) |>
    select(-data) |>
    unnest(performance) |>
    transmute(species_code = s,
              common_name = common_name,
              ppy_bin_midpoint,
              n_nonzero, directional_power, directional_error) |>
    bind_rows(ppy_binned)
  
  # abd bins, deciles by species
  p <- seq(0, 1, by = 0.1)
  p_midpoints <- (head(p, -1) + tail(p, -1)) / 2
  q <- quantile(t$abd, probs = p)
  idx <- cut(t$abd, q, include.lowest = TRUE) |>
    as.integer()
  t$abd_bin_midpoint <- p_midpoints[idx]
  abd_binned <- t |>
    # ppy magnitude
    mutate(ppy_magnitude = case_when(
      abs(simulated) <= 2 ~ "Weak",
      abs(simulated) <= 4 ~ "Moderate",
      abs(simulated) > 4 ~ "Strong",
    )) |>
    nest(data = !c(abd_bin_midpoint, ppy_magnitude)) |>
    mutate(performance = map(data, model_performance)) |>
    select(-data) |>
    unnest(performance) |>
    transmute(species_code = s,
              common_name = common_name,
              abd_bin_midpoint, ppy_magnitude,
              n_nonzero, directional_power, directional_error) |>
    bind_rows(abd_binned)
}


# power/error vs trend magnitude ----

plot_data <- ppy_binned |>
  filter(n_nonzero >= 25) |>
  select(common_name, ppy_bin_midpoint,
         directional_power, directional_error) |> 
  pivot_longer(cols = c(directional_power, directional_error)) |> 
  mutate(name = case_match(name,
                           "directional_power" ~ "Power",
                           "directional_error" ~ "Error"))
g <- ggplot(plot_data) +
  aes(x = ppy_bin_midpoint, y = value,
      color = common_name, linetype = name) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_brewer(palette = "Set1") +
  guides(linetype = guide_legend(nrow = 2),
         colour = guide_legend(nrow = 2)) +
  labs(x = "Trend magnitude [ |% change / year| ]",
       y = "Proportion of estimates",
       color = NULL, linetype = NULL) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
  theme_light() +
  theme(legend.position = "bottom")

ggsave(path(figures_dir, "figure-S09_power-error_ppy-bins.png"),
       g, width = 10, height = 8, units = "cm", scale = 2)
ggsave(path(figures_dir, "figure-S09_power-error_ppy-bins.tif"),
       g, width = 10, height = 8, units = "cm", scale = 2)


# power/error vs abundance ----

plot_data <- abd_binned |>
  filter(n_nonzero >= 25) |>
  select(common_name, abd_bin_midpoint, ppy_magnitude,
         directional_power, directional_error) |> 
  pivot_longer(cols = c(directional_power, directional_error)) |> 
  mutate(name = case_match(name,
                           "directional_power" ~ "Power",
                           "directional_error" ~ "Error"),
         ppy_magnitude = factor(ppy_magnitude,
                                levels = c("Weak", "Moderate", "Strong")))
g <- ggplot(plot_data) +
  aes(x = abd_bin_midpoint, y = value,
      color = common_name, shape = ppy_magnitude, linetype = name) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_brewer(palette = "Set1") +
  guides(linetype = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2),
         colour = guide_legend(nrow = 2)) +
  facet_wrap(~ common_name, nrow = 2) +
  labs(x = "Relative abundance quantile",
       y = "Proportion of estimates",
       shape = "Trend magnitude", color = NULL, linetype = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_light() +
  theme(legend.position = "bottom")

ggsave(path(figures_dir, "figure-S10_power-error_abd-bins.png"),
       g, width = 10, height = 8, units = "cm", scale = 2.3)
ggsave(path(figures_dir, "figure-S10_power-error_abd-bins.tif"),
       g, width = 10, height = 8, units = "cm", scale = 2.3)
