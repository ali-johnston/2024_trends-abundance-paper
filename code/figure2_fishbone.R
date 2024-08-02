library(arrow)
library(auk)
library(dplyr)
library(ebirdst)
library(fs)
library(ggplot2)
library(ggh4x)
library(glue)
library(readr)
library(rnaturalearth)
library(sf)
library(stringr)
library(tidyr)
library(purrr)

base_dir <- "/Users/Alison/Documents/REPOS/2024-trends-abundance-paper/"
data_dir <- path(base_dir, "data")
outputs_dir <- path(base_dir, "outputs")
figures_dir <- path(base_dir, "figures")

# read in trends estimates with breeding biomes and srd information
trends <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet()

# figure 2: min-max trends, all species ----

# trend quantiles across the whole range
trends_quantiles <- trends |>
  rename(abd_ppy = abd_ppy_median) |>
  group_by(species_code, season, breeding_biome) |>
  summarize(abd_ppy_rangewide = sum(abd * abd_ppy) / sum(abd),
            abd_ppy_min = min(abd_ppy),
            abd_ppy_5th = quantile(abd_ppy, 0.05) |> unname(),
            abd_ppy_25th = quantile(abd_ppy, 0.25) |> unname(),
            abd_ppy_median = quantile(abd_ppy, 0.50) |> unname(),
            abd_ppy_75th = quantile(abd_ppy, 0.75) |> unname(),
            abd_ppy_95th = quantile(abd_ppy, 0.95) |> unname(),
            abd_ppy_max = max(abd_ppy),
            .groups = "drop") |>
  arrange(abd_ppy_median) |>
  mutate(direction = ifelse(abd_ppy_rangewide > 0,
                            "increasing", "decreasing"),
         rank = row_number()) |>
  arrange(rank)


  fig_name <- "ebird-trends_min-max_iqr_median.png"
  trends_quantiles <- trends_quantiles |>
    mutate(plot_rank = rank)

fig_dir <- path(figures_dir, "figure2")
dir_create(fig_dir)
g <- ggplot(trends_quantiles) +
  aes(y = plot_rank) +
  # min-max decline
  geom_linerange(data = trends_quantiles |>
                   filter(abd_ppy_min < 0) |>
                   transmute(plot_rank, abd_ppy_min, abd_ppy_max = pmin(0, abd_ppy_max)),
                 aes(xmin = abd_ppy_min, xmax = abd_ppy_max),
                 color = "#FC9272") +
  # min-max increase
  geom_linerange(data = trends_quantiles |>
                   filter(abd_ppy_max > 0) |>
                   transmute(plot_rank, abd_ppy_min = pmax(0, abd_ppy_min), abd_ppy_max),
                 aes(xmin = abd_ppy_min, xmax = abd_ppy_max),
                 color = "#9ECAE1") +
  # iqr decline
  geom_linerange(data = trends_quantiles |>
                   filter(abd_ppy_25th < 0) |>
                   transmute(plot_rank, abd_ppy_25th, abd_ppy_75th = pmin(0, abd_ppy_75th)),
                 aes(xmin = abd_ppy_25th, xmax = abd_ppy_75th),
                 color = "#CB181D") +
  # iqr increase
  geom_linerange(data = trends_quantiles |>
                   filter(abd_ppy_75th > 0) |>
                   transmute(plot_rank, abd_ppy_25th = pmax(0, abd_ppy_25th), abd_ppy_75th),
                 aes(xmin = abd_ppy_25th, xmax = abd_ppy_75th),
                 color = "#2171B5", linewidth = 0.5) +
  # rangewide trend
  geom_point(aes(x = abd_ppy_median), size = 0.1) +
  scale_x_continuous(breaks = seq(-20, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(xlim = c(-15, 15), ylim = c(-5, 496)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
        axis.text = element_text(size = rel(1.7)),
        axis.text.y = element_blank())
ggsave(path(fig_dir, fig_name), g,
       width = 8, height = 12, bg = "white")


