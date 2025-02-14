library(arrow)
library(dplyr)
library(fs)
library(ggplot2)
library(readr)
library(tidyr)

data_dir <- "data"
figures_dir <- path("figures", "figure-S11")
dir_create(figures_dir)


# load ----

# performance metrics binned by trend magnitude
ppy_binned <- path(data_dir,
                   "ebird-trends_ppy-binned-performance_2021.parquet") |>
  read_parquet() |>
  select(species_code, season,
         ppy_bin_midpoint,
         n_nonzero,
         directional_power = directional_power,
         directional_error = directional_error_nz) |> 
  filter(n_nonzero >= 25) |>
  drop_na()
# performance metrics binned by abundance
abd_binned <- path(data_dir,
                   "ebird-trends_abd-binned-performance_2021.parquet") |>
  read_parquet() |>
  select(species_code, season,
         abd_bin_midpoint, n_nonzero,
         directional_power = directional_power,
         directional_error = directional_error_nz) |> 
  filter(n_nonzero >= 25) |>
  drop_na()


# figure S11a: power vs ppy ----

# quartiles of power
ptiles <- ppy_binned |>
  group_by(ppy_bin_midpoint) |>
  summarize(p25 = quantile(directional_power, 0.25),
            p50 = quantile(directional_power, 0.5),
            p75 = quantile(directional_power, 0.75),
            .groups = "drop") |>
  pivot_longer(c(p25, p50, p75),
               names_to = "percentile",
               values_to = "directional_power")

# plot
gga <- ggplot(ppy_binned) +
  aes(x = ppy_bin_midpoint, y = directional_power) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(aes(group = species_code),
            color = scales::alpha("black", 0.03)) +
  geom_line(data = ptiles, mapping = aes(group = percentile),
            color = "red") +
  geom_vline(xintercept = 3.3, linetype = "dashed") +
  geom_vline(xintercept = 6.7, linetype = "dashed") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Trend magnitude [ |% change / year| ]",
       y = "Power", color = NULL) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 15)) +
  theme_light() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.75)))
ggsave(path(figures_dir, "figure-S11a_power-vs-ppy.png"),
       plot = gga, width = 10, height = 8, scale = 0.6)
ggsave(path(figures_dir, "figure-S11a_power-vs-ppy.tif"),
       plot = gga, width = 10, height = 8, scale = 0.6, dpi = 600)


# figure S11b: power vs abd ----

# quartiles of power
ptiles <- abd_binned |>
  group_by(abd_bin_midpoint) |>
  summarize(p25 = quantile(directional_power, 0.25),
            p50 = quantile(directional_power, 0.5),
            p75 = quantile(directional_power, 0.75),
            .groups = "drop") |>
  pivot_longer(c(p25, p50, p75),
               names_to = "percentile",
               values_to = "directional_power")

# plot
ggb <- ggplot(abd_binned) +
  aes(x = abd_bin_midpoint, y = directional_power) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(aes(group = species_code),
            color = scales::alpha("black", 0.03)) +
  geom_line(data = ptiles, mapping = aes(group = percentile),
            color = "red") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Relative abundance quantile",
       y = "Power", color = NULL) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  theme_light() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.75)))
ggsave(path(figures_dir, "figure-S11b_power-vs-abd.png"),
       plot = ggb, width = 10, height = 8, scale = 0.6)
ggsave(path(figures_dir, "figure-S11b_power-vs-abd.tif"),
       plot = ggb, width = 10, height = 8, scale = 0.6, dpi = 600)


# figure S11c: error vs ppy ----

# quartiles of error
ptiles <- ppy_binned |>
  group_by(ppy_bin_midpoint) |>
  summarize(p25 = quantile(directional_error, 0.25),
            p50 = quantile(directional_error, 0.5),
            p75 = quantile(directional_error, 0.75),
            .groups = "drop") |>
  pivot_longer(c(p25, p50, p75),
               names_to = "percentile",
               values_to = "directional_error")

# plot
ggc <- ggplot(ppy_binned) +
  aes(x = ppy_bin_midpoint, y = directional_error) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(aes(group = species_code),
            color = scales::alpha("black", 0.03)) +
  geom_line(data = ptiles, mapping = aes(group = percentile),
            color = "red") +
  geom_vline(xintercept = 3.3, linetype = "dashed") +
  geom_vline(xintercept = 6.7, linetype = "dashed") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Trend magnitude [ |% change / year| ]",
       y = "Error", color = NULL) +
  coord_cartesian(ylim = c(0, 0.6), xlim = c(0, 15)) +
  theme_light() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.75)))
ggsave(path(figures_dir, "figure-S11c_error-vs-ppy.png"),
       plot = ggc, width = 10, height = 8, scale = 0.6)
ggsave(path(figures_dir, "figure-S11c_error-vs-ppy.tif"),
       plot = ggc, width = 10, height = 8, scale = 0.6, dpi = 600)


# figure S11c: error vs ppy ----

# quartiles of error
ptiles <- abd_binned |>
  group_by(abd_bin_midpoint) |>
  summarize(p25 = quantile(directional_error, 0.25),
            p50 = quantile(directional_error, 0.5),
            p75 = quantile(directional_error, 0.75),
            .groups = "drop") |>
  pivot_longer(c(p25, p50, p75),
               names_to = "percentile",
               values_to = "directional_error")

# plot
ggd <- ggplot(abd_binned) +
  aes(x = abd_bin_midpoint, y = directional_error) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(aes(group = species_code),
            color = scales::alpha("black", 0.03)) +
  geom_line(data = ptiles, mapping = aes(group = percentile),
            color = "red") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Relative abundance quantile",
       y = "Error", color = NULL) +
  coord_cartesian(ylim = c(0, 0.6), xlim = c(0, 1)) +
  theme_light() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.75)))
ggsave(path(figures_dir, "figure-S11d_error-vs-abd.png"),
       plot = ggd, width = 10, height = 8, scale = 0.6)
ggsave(path(figures_dir, "figure-S11d_error-vs-abd.tif"),
       plot = ggd, width = 10, height = 8, scale = 0.6, dpi = 600)
