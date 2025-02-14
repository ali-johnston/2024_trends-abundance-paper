library(arrow)
library(dplyr)
library(fs)
library(ggplot2)
library(ggh4x)
library(tidyr)

data_dir <- "data"
outputs_dir <- "outputs"
figures_dir <- path("figures", "figure-S03-S04")
dir_create(figures_dir)

# load data ----

# range-wide trends for 100 folds
trends_rw_folds <- path(data_dir,
                        "ebird-trends_range-wide_folds_2021.parquet") |>
  read_parquet()

# order biomes by amount of decline
biome_order <- trends_rw_folds |>
  group_by(breeding_biome) |>
  summarise(prop_decline_range = median(prop_decline_range), .groups ="drop") |>
  arrange(desc(prop_decline_range)) |>
  pull(breeding_biome)
trends_rw_folds$breeding_biome <- factor(trends_rw_folds$breeding_biome,
                                         levels = biome_order,
                                         ordered = TRUE)

# reformat to show abd_ppy and prop_decline_range on same plot
trends_rw_folds <- trends_rw_folds |>
  select(breeding_biome, species_code, season, fold,
         abd_ppy, prop_decline_range) |>
  pivot_longer(cols = c(abd_ppy, prop_decline_range),
               names_to = "type", values_to = "metric") |>
  mutate(type_label = if_else(type == "abd_ppy",
                              "Range-wide population trend",
                              "Proportion of range with declines"))

# figS3: density plot ----

# create a fake little data frame to allow different vertical lines
# in each column:
vline_df <- trends_rw_folds |>
  distinct(type, type_label, breeding_biome) |>
  mutate(x_int = ifelse(type == "abd_ppy", 0, 0.5))

# density plots by biome and fold
g <- ggplot(trends_rw_folds, aes(x = metric)) +
  geom_vline(data = vline_df, linewidth = 0.25,
             aes(xintercept = x_int),
             colour = alpha("black", 0.4)) +
  geom_density(aes(group = as.factor(fold)),
               color = alpha("firebrick", 0.05),
               show.legend = FALSE) +
  facet_grid2(vars(breeding_biome), vars(type_label),
              scales = "free", independent = "y",
              switch = "both") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        strip.background =element_rect(colour = NA, fill="white"),
        strip.placement = "outside") +
  ylab("Relative number of species") +
  xlab(NULL)
ggsave(path(figures_dir, "figure-S03_proportion-cells_density-plot.png"), g,
       width = 7, height = 8, bg = "white")
ggsave(path(figures_dir, "figure-S03_proportion-cells_density-plot.tif"), g,
       width = 7, height = 8, bg = "white")


# figS4: density plot by season ----

g <- ggplot(trends_rw_folds |>
              mutate(season = ifelse(season == "resident", "breeding", season)),
            aes(x = metric)) +
  geom_vline(data = vline_df, linewidth = 0.25,
             aes(xintercept = x_int),
             colour = alpha("black", 0.4)) +
  geom_density(aes(group = interaction(fold, season), color = season),
               show.legend = FALSE) +
  scale_color_manual(values = c(alpha("firebrick", 0.05),
                                alpha("steelblue", 0.05)),
                     breaks = c("breeding", "nonbreeding")) +
  facet_grid2(vars(breeding_biome), vars(type_label),
              scales = "free", independent = "y",
              switch = "both") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        strip.background =element_rect(colour = NA, fill="white"),
        strip.placement = "outside") +
  ylab("Relative number of species") +
  xlab(NULL)
ggsave(path(figures_dir, "figure-S04_proportion-cells_density-plot_by-season.png"), g,
       width = 7, height = 8, bg = "white")
ggsave(path(figures_dir, "figure-S04_proportion-cells_density-plot_by-season.tif"), g,
       width = 7, height = 8, bg = "white")
