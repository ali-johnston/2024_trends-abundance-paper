#########################################################
## fit the linear mixed model for abundance and 
## create figure 3

library(tidyverse)
library(lme4)
library(arrow)
library(fs)

base_dir <- "/Users/Alison/Documents/REPOS/2024-trends-abundance-paper/"
data_dir <- path(base_dir, "data")
outputs_dir <- path(base_dir, "outputs")
figures_dir <- path(base_dir, "figures/figure3")
dir_create(figures_dir)

#########################################################
## load data

# # species lookup
species <- read_csv("data/master_species_list_495.csv", na = "") |>
  select(species_code, breeding_biome)

# trends estimates
trends <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet() |>
  mutate(log_abd = log10(abd),
#         log_distance_to_edge_km = log10(distance_to_edge_km),
         breeding_biome = factor(breeding_biome),
         species_code = factor(species_code),
         weight = 1 / abd_ppy_var)

# minimum and maximum abundance for each species and each biome
abd_range <- trends |>
  dplyr::select(species_code, breeding_biome, abd) |>
  group_by(species_code, breeding_biome) |>
  summarise(max_abd = max(abd), min_abd = min(abd), .groups = "drop") |>
  mutate(min_log_abd = log10(min_abd), max_log_abd = log10(max_abd))

abd_range_biome <- abd_range |>
  dplyr::select(breeding_biome, min_log_abd, max_log_abd) |>
  group_by(breeding_biome) |>
  summarise(min_log_abd = min(min_log_abd),
            max_log_abd = max(max_log_abd), .groups = "drop")


#########################################################
## Mixed model with abundance and 1/var as weight: 

# mixed model
m1 <- lmer(abd_ppy_median ~ log_abd + (1 + log_abd|breeding_biome/species_code),
           weights = weight,
           data = trends)

mod <- m1
mod_tag <- "abd"
var_to_plot <- "log_abd" # "log_distance_to_edge_km"

#########################################################
## extract slope and intercept estimates and save

re <- ranef(mod)
fi <- fixef(mod)

biome_coef_0 <- re[["breeding_biome"]] |>
  rownames_to_column(var = "breeding_biome") |>
  rename(biome_int = "(Intercept)",
         biome_log_abd = "log_abd")

spec_coef <- re[["species_code:breeding_biome"]] |>
  rownames_to_column(var = "spec_biome") |>
  rename(spec_int = "(Intercept)",
         spec_log_abd = "log_abd") |>
  separate_wider_delim(spec_biome, delim = ":",
                       names = c("species_code", "breeding_biome")) |>
  inner_join(biome_coef_0, by = join_by(breeding_biome)) |>
  mutate(overall_int = fi[1] + spec_int + biome_int) |>
  mutate(overall_log_abd = fi[2] + spec_log_abd + biome_log_abd) |>
  inner_join(abd_range, by = join_by(species_code, breeding_biome)) |>
  arrange(overall_log_abd)

biome_coef <- biome_coef_0 |>
  mutate(overall_int = fi[1] + biome_int,
         overall_log_abd = fi[2] + biome_log_abd) |>
  inner_join(abd_range_biome, by = join_by(breeding_biome)) |>
  arrange(overall_log_abd)

write_csv(spec_coef, path(outputs_dir, paste0("mixed-model_species-coefs_", mod_tag, ".csv")))
write_csv(biome_coef, path(outputs_dir, paste0("mixed-model_biome-coefs_", mod_tag, ".csv")))



#########################################################
##  figure 3a: histogram of slopes ----

# max slope for limits
max_slope <- ceiling(max(abs(spec_coef$overall_log_abd))*10)/10

ggplot(spec_coef, aes(x = overall_log_abd)) +
  geom_histogram(breaks = seq(-1*max_slope, max_slope, by = 0.1),
                 fill = alpha("steelblue3", 0.8),
                 colour = "white", show.legend = FALSE) +
  #scale_fill_manual(values = c(alpha("steelblue3", 0.8), alpha("steelblue3", 0.4))) +
  theme_bw() +
  geom_vline(xintercept = 0) +
  ylab("Number of species") +
  xlab("Mixed model slope") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_loc <- path(figures_dir, paste0("fig3a_mixed-model_slope-histogram_", mod_tag, ".png"))
ggsave(plot_loc, width = 14, height = 10, units = "cm")



#########################################################
## figure 3b: trend vs. abundance relationship from random effects ----

plot_1_line <- function(data, col = alpha("black", 0.1), ...){
  x <- seq(0, 1, by = 0.1)
  x2 <- seq(data$min_log_abd[1], data$max_log_abd[1], length.out = 11)
  lines(x, x2*data$overall_log_abd[1] + data$overall_int[1], col = col, ...)
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

plot_loc <- path(figure_dir, paste0("fig3b_mixed-model_slopes_", mod_tag, ".png"))
png(plot_loc, width = 14, height = 10, units = "cm", pointsize = 9, res = 600)

par(mfrow=c(2, 3), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,5,2,0))


for(i in 1:length(biomes)) {
  ylabel <- ifelse(i %in% c(1, 4), "Population trend", "")
  plot(0, 0, col = "white", ylim = c(-12, 12), xlim = c(0, 1),
       main = "", xaxt = "n", xlab = "", yaxt = "n", ylab = ylabel,
       xpd = NA, cex.lab = 1.5, las = 1)
  if(i %in% c(1,4)) axis(side = 2, at = c(-10, 0, 10), cex.axis = 1.5, las = 1)
  if(i > 3) axis(side = 1, at = c(0, 1), labels = c("0", "max"), cex.axis = 1.5, las = 1)
  abline(h = 0, col = alpha("firebrick", 0.42), lwd = 2)

  # lines for each species in the biome
  plot_multi(filter(spec_coef, breeding_biome == biomes[i]))

  if(i == 5) text(x = 0.5, y = -20,
                  labels = "Relative abundance within species", xpd = NA, cex = 1.5)
  text(x = 0, y = 11, pos = 4, labels = biome_names_tidy[i], cex = 1.5)

  # add average biome line
  biome_avg <- filter(biome_coef, breeding_biome == biomes[i])
  plot_1_line(biome_avg, col = "black", lwd = 2)
}

dev.off()



#########################################################
## individual species plots

focal_species <- c("rufhum", "balori")

species_names <- path(data_dir, "master_species_list_495.csv") |>
  read_csv() |>
  filter(species_code %in% focal_species)

for(i in 1:nrow(species_names)){

  # collect species data for plots
  this_species_code <- species_names$species_code[i]
  this_species_name <- gsub("'", "", species_names$common_name[i], fixed = TRUE)
  this_species_name <- gsub(" ", "", this_species_name, fixed = TRUE)
  this_species_coef <- spec_coef |> filter(species_code == this_species_code)
  this_species_trends <- trends |> filter(species_code == this_species_code)

  if(nrow(this_species_coef)==0) {
    # print out the names of any species without data. 
    # if all is working well, there should be no species names printed 
    print(paste("species", this_species_code, this_species_name))
  } else {

  # set up plot (to hopefully make pretty for many species)
  xlimits <- range(this_species_trends$log_abd)
  la_vals <- seq(-5, 5, by = 1)

  if(xlimits[2] - xlimits[1])
  la_vals <- seq(-6, 6, by = 2)

  ylimits <- range(this_species_trends$abd_ppy_median)
  if(ylimits[2]< -10) ylimits[2] <- -10
  if(ylimits[1]> 10) ylimits[1] <- 10
  ylimits <- c(-10, 10)


  formatted_slope <- paste0(ifelse(this_species_coef$overall_log_abd<0, "n", "p"), "_",
           format(round(abs(this_species_coef$overall_log_abd), digits = 3), nsmall = 3))
  species_plot_name <- paste0(formatted_slope, "_", this_species_code, "_", this_species_name)

  plot_loc <- path(figures_dir, paste0("species_scatter_", species_plot_name, ".png"))
  png(plot_loc, width = 8, height = 8, units = "cm", pointsize = 9, res = 600)
  par(mar = c(5, 5, 3, 3))

  # create plot
  plot(0, 0, col = "white", 
    xlim = xlimits, 
    xaxt = "n", xlab = "Relative abundance",
    ylim = ylimits,
    yaxt = "n", ylab = "Population trend", 
    main = this_species_name, 
    mgp = c(1.5, 0.7, 0))
  axis(side = 1, at = la_vals, labels = 10^la_vals, tcl = -0.3, mgp = c(1.5, 0.5, 0))
  axis(side = 2, at = c(-10, 0, 10), las = 1, tcl = -0.3, mgp = c(1.5, 0.7, 0))

  abline(h = 0, col = alpha("firebrick", 0.2), lwd = 2)

  points(this_species_trends$log_abd, this_species_trends$abd_ppy_median, 
        pch = 16, col = alpha("steelblue", 0.3),
        cex = sqrt(this_species_trends$weight))

  abline(this_species_coef$overall_int, this_species_coef$overall_log_abd, col = "steelblue", lwd= 1.5)

  dev.off()
}

}



