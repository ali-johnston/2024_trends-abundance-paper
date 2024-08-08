library(arrow)
library(dplyr)
library(fs)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyr)
library(auk)

base_dir <- "/Users/Alison/Documents/REPOS/2024_trends-abundance-paper/"
data_dir <- path(base_dir, "data")
figures_dir <- path(base_dir, "figures", "figure_S9-S10")
dir_create(figures_dir)

tax <- select(auk::ebird_taxonomy, species_code, common_name)

# trends data
trends <- path(data_dir,
               "ebird-trends_simulations_focal-species_2021.parquet") |>
  read_parquet() |>
  mutate(direction_sim = sign(simulated),
         direction_est = sign(estimated),
         correct_direction = direction_sim == direction_est,
         correct_direction_nz = nonzero & correct_direction)

# attach abundnace
trends <- path(data_dir,
            "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet() |>
  select(species_code, srd_id, abd) |>
  inner_join(trends, by = join_by(species_code, srd_id))

# power vs trend ----
bin_trends <- function(data, bin_width = 0.5, bin_column = "estimated") {
  stopifnot(is.data.frame(data), bin_column %in% names(data))
  stopifnot(is.numeric(bin_width), length(bin_width) == 1, bin_width > 0)

  # assign bin index
  data$bin_index <- as.integer(abs(data[[bin_column]]) %/% bin_width)

  # create bin midpoint lookup table
  midpoints <- data.frame(bin_index = seq(0, max(data$bin_index), by = 1))
  midpoints$bin_midpoint <- bin_width * midpoints$bin_index + bin_width / 2
  attr(data, "midpoints") <- midpoints

  return(data)
}

binned_trends <- trends |>
  bin_trends() |>
  mutate(trends_bin_index = bin_index)
trends_midpoints <- attr(binned_trends, "midpoints") |>
  select(trends_bin_index = bin_index,
         trends_bin_midpoint = bin_midpoint)
power_trends <- binned_trends |>
  group_by(species_code, season, trends_bin_index) |>
  summarise(directional_power = mean(correct_direction_nz),
            .groups = "drop") |>
  inner_join(tax, by = join_by(species_code)) |>
  inner_join(trends_midpoints, by = join_by(trends_bin_index)) |>
  drop_na() |>
  select(species_code, season, common_name,
         trends_bin_midpoint, directional_power)

gs8a <- ggplot(power_trends) +
  aes(x = trends_bin_midpoint, y = directional_power,
      group = common_name, color = common_name) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                     labels = scales::label_percent()) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Trend Magnitude [ | % change / year | ]",
       y = "Power", color = NULL) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 10)) +
  theme_light() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.75)))
ggsave(path(figures_dir, "figure_s9_power_trend-magnitude.png"), plot = gs8a,
       width = 10, height = 8, scale = 0.6)

ggsave(path(figures_dir, "figure_s9_power_trend-magnitude.tiff"), plot = gs8a,
       width = 10, height = 8, scale = 0.6)

# power vs. abundance----
binned_abd <- NULL
for (s in unique(trends$species_code)) {
  q <- seq(0, 1, by = 0.1)
  trends_sp <- filter(trends, species_code == s)
  abd_bins <- quantile(trends_sp$abd, q)
  midpoints <- cbind(q[-length(q)], q[-1]) |>
    rowMeans() |>
    unname()
  trends_sp$abd_quantile <- cut(trends_sp$abd,
                                breaks = abd_bins,
                                labels = midpoints)
  binned_abd <- trends_sp |>
    mutate(abd_bin_index = as.integer(abd_quantile),
           abd_quantile = as.numeric(as.character(abd_quantile))) |>
    bind_rows(binned_abd)
}
# add trend magnitude bin
binned_abd <- binned_abd |>
  mutate(trend_bin = case_when(
    abs(estimated) <= 1 ~ "Weak",
    abs(estimated) >= 2 & abs(estimated) <= 4 ~ "Moderate",
    abs(estimated) >= 5 ~ "Strong",
    .default = NA_character_),
    trend_bin = factor(trend_bin, levels = c("Weak", "Moderate", "Strong")))
power_abd <- binned_abd |>
  group_by(species_code, season, abd_quantile, abd_bin_index) |>
  summarise(directional_power = mean(correct_direction_nz),
            .groups = "drop") |>
  inner_join(tax, by = join_by(species_code)) |>
  select(species_code, season, common_name,
         abd_quantile, directional_power) |>
  drop_na()

gs8b <- ggplot(power_abd) +
  aes(x = abd_quantile, y = directional_power,
      group = common_name, color = common_name) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                     labels = scales::label_percent()) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                     labels = scales::label_percent()) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Relative Abundance Quantile",
       y = "Power", color = NULL) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  theme_light() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.75)))
ggsave(path(fig_dir, "s9b_power_abundance.png"), plot = gs8b,
       width = 10, height = 8, scale = 0.6)
# stack figures
ggsave(path(fig_dir, "s9_power_stacked.tif"),
       plot = gs8a / gs8b,
       width = 10, height = 16, scale = 0.6, dpi = 600)


# power vs abundance group by trend magnitude ----
power_abd_trend <- binned_abd |>
  group_by(species_code, season, abd_quantile, trend_bin) |>
  summarise(directional_power = mean(correct_direction_nz),
            .groups = "drop") |>
  inner_join(tax, by = join_by(species_code)) |>
  select(species_code, season, common_name,
         abd_quantile, trend_bin, directional_power) |>
  drop_na()
gs9 <- ggplot(power_abd_trend) +
  aes(x = abd_quantile, y = directional_power,
      color = common_name, shape = trend_bin) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.20),
                     labels = scales::label_percent()) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = scales::label_percent()) +
  scale_color_brewer(palette = "Set1", guide = "none") +
  labs(x = "Relative Abundance Quantile",
       y = "Power", color = NULL,
       shape = "Trend Magnitude") +
  facet_wrap(~ common_name) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  theme_light() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.75)))
ggsave(path(figures_dir, "s10_power_interaction.png"), plot = gs9,
       width = 10, height = 10, scale = 0.6)
ggsave(path(figures_dir, "s10_power_interaction.tiff"), plot = gs9,
       width = 10, height = 10, scale = 0.6, dpi = 600)

