library(arrow)
library(dplyr)
library(fs)
library(ggplot2)
library(readr)
library(tidyr)

data_dir <- "data"
figures_dir <- path("figures", "figure-S12")
dir_create(figures_dir)


# load ----

# performance metrics binned by trend magnitude
ppy_binned <- path(data_dir,
                   "ebird-trends_ppy-binned-performance_2021.parquet") |>
  read_parquet() |>
  group_by(ppy_bin_midpoint) |> 
  summarize(n = sum(n_nonzero / prop_nonzero, na.rm = TRUE),
            n_correct = sum(n_nonzero * directional_power_nz, na.rm = TRUE),
            n_incorrect = sum(n_nonzero * directional_error_nz, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(Correct = n_correct / n,
         Incorrect = n_incorrect / n) |> 
  pivot_longer(cols = c(Correct, Incorrect))


# figure S12: power-error histogram ----

g <- ggplot(ppy_binned) +
  aes(x = ppy_bin_midpoint, y = value, fill = name) +
  geom_bar(position="stack", stat="identity") +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Trend magnitude [ |% change / year| ]",
       y = "Percent of all estimates",
       fill = "Direction of non-zero estimates") +
  coord_cartesian() +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(path(figures_dir, "figure-S12_power-error-histogram.png"),
       g, width = 10, height = 8, scale = 0.6)
ggsave(path(figures_dir, "figure-S12_power-error-histogram.tif"),
       g, width = 10, height = 8, scale = 0.6, dpi = 600)
