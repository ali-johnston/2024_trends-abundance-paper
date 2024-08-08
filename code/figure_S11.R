library(arrow)
library(dplyr)
library(fs)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyr)


base_dir <- "/Users/Alison/Documents/REPOS/2024_trends-abundance-paper/"
data_dir <- path(base_dir, "data")
figures_dir <- path(base_dir, "figures", "figure_S11")
dir_create(figures_dir)

# all species ----
binned_perf <- path(data_dir,
                    "ebird-trends_binned-performance_2021.parquet") |>
  read_parquet() |>
  filter(n_nonzero >= 25) |>
  drop_na()
ptiles <- binned_perf |>
  group_by(abd_ppy_midpoint) |>
  summarize(p25 = quantile(directional_power, 0.25),
            p50 = quantile(directional_power, 0.5),
            p75 = quantile(directional_power, 0.75),
            .groups = "drop") |>
  pivot_longer(c(p25, p50, p75),
               names_to = "percentile",
               values_to = "directional_power")
gs10 <- ggplot(binned_perf) +
  aes(x = abd_ppy_midpoint, y = directional_power) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(aes(group = species_code),
            color = scales::alpha("black", 0.03)) +
  geom_line(data = ptiles, mapping = aes(group = percentile),
            color = "red") +
  geom_vline(xintercept = 3.3, linetype = "dashed") +
  geom_vline(xintercept = 6.7, linetype = "dashed") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Trend Magnitude [ | % change / year | ]",
       y = "Power", color = NULL) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 15)) +
  theme_light() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.75)))
ggsave(path(figures_dir, "fig_S11_power_trend-magnitude_all-species.png"),
       plot = gs10, width = 10, height = 8, scale = 0.6)
ggsave(path(figures_dir, "fig_S11_power_trend-magnitude_all-species.tiff"),
       plot = gs10, width = 10, height = 8, scale = 0.6, dpi = 600)
