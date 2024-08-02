base_dir <- "/Users/Alison/Documents/REPOS/2024-trends-abundance-paper/"
data_dir <- path(base_dir, "data")
outputs_dir <- path(base_dir, "outputs")
figures_dir <- path(base_dir, "figures/figure_S11_S12")
dir.create(figures_dir)

library(tidyverse)
library(lme4)
library(arrow)
library(fs)


#########################################################
## load data

# # species lookup
species <- read_csv("data/master_species_list_495.csv", na = "") |>
  select(species_code, breeding_biome)

# trends estimates
trends <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet() |>
  mutate(log_abd = log10(abd),
         log_distance_to_edge_km = log10(distance_to_edge_km),
         breeding_biome = factor(breeding_biome),
         species_code = factor(species_code),
         weight = 1 / abd_ppy_var)


# minimum and maximum distance to edge for each species and each biome
# used for predictions and plots

dist_range <- trends |>
  dplyr::select(species_code, breeding_biome, distance_to_edge_km) |>
  group_by(species_code, breeding_biome) |>
  summarise(max_dist = max(distance_to_edge_km), min_dist = min(distance_to_edge_km), .groups = "drop") |>
  mutate(min_log_dist = log10(min_dist), max_log_dist = log10(max_dist))

dist_range_biome <- dist_range |>
  dplyr::select(breeding_biome, min_log_dist, max_log_dist) |>
  group_by(breeding_biome) |>
  summarise(min_log_dist = min(min_log_dist),
            max_log_dist = max(max_log_dist), .groups = "drop")


#########################################################
## Mixed model with distance to edge and 1/var as weight: 

# mixed model
m3 <- lmer(abd_ppy_median ~ log_distance_to_edge_km + (1 + log_distance_to_edge_km|breeding_biome/species_code),
           weights = weight,
           data = trends)

mod <- m3
mod_tag <- "edge"
x_axis <- "log_distance_to_edge_km"

#########################################################
## extract slope and intercept estimates and save

re <- ranef(mod)
fi <- fixef(mod)

biome_coef_0 <- re[["breeding_biome"]] |>
  rownames_to_column(var = "breeding_biome") |>
  rename(biome_int = "(Intercept)",
         biome_log_distance_to_edge_km = "log_distance_to_edge_km")

spec_coef <- re[["species_code:breeding_biome"]] |>
  rownames_to_column(var = "spec_biome") |>
  rename(spec_int = "(Intercept)",
         spec_log_distance_to_edge_km = "log_distance_to_edge_km") |>
  separate_wider_delim(spec_biome, delim = ":",
                       names = c("species_code", "breeding_biome")) |>
  inner_join(biome_coef_0, by = join_by(breeding_biome)) |>
  mutate(overall_int = fi[1] + spec_int + biome_int) |>
  mutate(overall_log_distance_to_edge_km = fi[2] + spec_log_distance_to_edge_km + biome_log_distance_to_edge_km) |>
  inner_join(dist_range, by = join_by(species_code, breeding_biome)) |>
  arrange(overall_log_distance_to_edge_km)

biome_coef <- biome_coef_0 |>
  mutate(overall_int = fi[1] + biome_int,
         overall_log_distance_to_edge_km = fi[2] + biome_log_distance_to_edge_km) |>
  inner_join(dist_range_biome, by = join_by(breeding_biome)) |>
  arrange(overall_log_distance_to_edge_km)

write_csv(spec_coef, path(outputs_dir, paste0("mixed-model_species-coefs_", mod_tag, ".csv")))
write_csv(biome_coef, path(outputs_dir, paste0("mixed-model_biome-coefs_", mod_tag, ".csv")))



#########################################################
##  figure 2a: histogram of slopes ----

# max slope for limits
max_slope <- ceiling(max(abs(spec_coef$overall_log_distance_to_edge_km))*10)/10

ggplot(spec_coef, aes(x = overall_log_distance_to_edge_km)) +
  geom_histogram(breaks = seq(-1*max_slope, max_slope, length.out =),
                 fill = alpha("steelblue3", 0.8),
                 colour = "white", show.legend = FALSE) +
  #scale_fill_manual(values = c(alpha("steelblue3", 0.8), alpha("steelblue3", 0.4))) +
  theme_bw() +
  geom_vline(xintercept = 0) +
  ylab("Number of species") +
  xlab("Mixed model slope") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot_loc <- path(figures_dir, paste0("fig11_mixed-model_slope-histogram_", mod_tag, ".png"))
ggsave(plot_loc, width = 14, height = 10, units = "cm")



#########################################################
## trend vs. distance to edge relationships from random effects ----

plot_1_line <- function(data, col = alpha("black", 0.1), ...){
  x <- seq(0, 1, by = 0.1)
  x2 <- seq(data$min_log_dist[1], data$max_log_dist[1], length.out = 11)
  lines(x, x2*data$overall_log_distance_to_edge_km[1] + data$overall_int[1], col = col, ...)
}
plot_multi <- function(data, col = alpha("black", 0.1)){
  for(i in 1:nrow(data)){
    plot_1_line(data[i,], col = col)
  }
}

# change breeding biome names to match paper
biomes <- names(table(spec_coef$breeding_biome))
biome_names_tidy <- c("Arctic tundra", "Aridland", "Forest",
      "Grassland", "Habitat generalists", "Wetland & Coast")


plot_loc <- path(figures_dir, paste0("fig12_mixed-model_slopes_", mod_tag, ".png"))
png(plot_loc, width = 14, height = 10, units = "cm", pointsize = 9, res = 600)
par(mfrow=c(2, 3), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,5,2,0))
biomes <- names(table(spec_coef$breeding_biome))

for (i in seq_along(biomes)) {
  ylabel <- ifelse(i %in% c(1, 4), "Population trend", "")
  plot(0, 0, col = "white", ylim = c(-12, 12), xlim = c(0, 1),
       main = "", xaxt = "n", xlab = "", yaxt = "n", ylab = ylabel,
       xpd = NA, cex.lab = 1.5)
  if(i %in% c(1,4)) axis(side = 2, at = c(-10, -5, 0, 5, 10), cex.axis = 1.5)
  if(i > 3) axis(side = 1, at = c(0, 1), labels = c("0", "max"), cex.axis = 1.5)
  abline(h = 0, col = alpha("firebrick", 0.42), lwd = 2)

  # lines for each species in the biome
  plot_multi(filter(spec_coef, breeding_biome == biomes[i]))

  if(i == 5) text(x = 0.5, y = -17,
                  labels = "Relative distance to edge within species", xpd = NA, cex = 1.5)
  text(x = 0, y = 11, pos = 4, labels = biome_names_tidy[i], cex = 1.5)

  # add average biome line
  biome_avg <- filter(biome_coef, breeding_biome == biomes[i])
  plot_1_line(biome_avg, col = "black", lwd = 2)
}

dev.off()


