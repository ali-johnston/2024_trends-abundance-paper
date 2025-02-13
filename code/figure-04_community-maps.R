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

data_dir <- "data"
outputs_dir <- "outputs"
figures_dir <- path("figures", "figure-04")
dir_create(figures_dir)

# read in trends estimates with breeding biomes and srd information
trends <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet()

# srd <- path(data_dir, "srd_27km_year.parquet") |>
srd <- path(data_dir, "srd_27km.parquet") |>
  read_parquet()

# map plotting set-up
# focal region
crs <- st_crs("+proj=laea +lon_0=-100 +lat_0=50")
basemap <- list()
basemap_gpkg <- path(data_dir, "basemap.gpkg")
for (n in st_layers(basemap_gpkg)$name) {
  basemap[[n]] <- read_sf(basemap_gpkg, layer = n) |>
    st_transform(crs) |>
    st_geometry()
}

# lower 48
region <- ne_states(iso_a2 = c("US", "CA"), returnclass = "sf") |>
  filter(!postal %in% c("AK", "HI"), iso_a2 %in% c("US") | postal == "NS") |>
  st_transform(crs)
bb <- st_bbox(c(xmin = -2040211, ymin = -3656994,
                xmax = 3060613, ymax = 1443830)) |>
  st_as_sfc()

# theme
max_radius <- 0.95 * 26665.26 / 2
bg_col <- "#ffffff"
land_col <- "#ffffff"
land_border_col <- "#888888"
line_col <- "#888888"
line_width <- 2.5

# biome level trends
trends_biomes <- trends |>
  filter(longitude > -165) |>
  group_by(breeding_biome, breeding_biome_label, srd_id) |>
  summarise(abd_ppy_mean = mean(abd_ppy_median),
            decline_prop = mean(abd_ppy_median < 0),
            n_species = n(),
            .groups = "drop") |>
  inner_join(srd, by = "srd_id") |>
  select(breeding_biome, breeding_biome_label, srd_id, longitude, latitude,
         abd_ppy_mean, decline_prop, n_species)

trends_biomes_sf <- trends_biomes |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs = crs)

# bounding box for maps
bb <- st_bbox(c(xmin = -2040211, ymin = -3656994,
                xmax = 3060613, ymax = 1443830),
              crs = crs) |>
  st_as_sfc()

# ppy bins
ppy_bins <- seq(0, 4, by = 0.5)
ppy_bins[1] <- min(abs(trends_biomes_sf$abd_ppy_mean)) / 10
ppy_bins[length(ppy_bins)] <- max(abs(trends_biomes_sf$abd_ppy_mean)) * 10
ppy_bins <- c(-rev(sort(ppy_bins)), sort(ppy_bins))
pal <- ebirdst_palettes(length(ppy_bins) - 1, type = "trends")
cols <- cut(trends_biomes_sf$abd_ppy_mean, ppy_bins)
trends_biomes_sf$abd_ppy_cols <- pal[match(cols, levels(cols))]

# proportion of declines bins
decline_bins <- seq(0, 1, by = 0.05)
decline_bins[1] <- -Inf
decline_bins[length(decline_bins)] <- Inf
# assign colors based on proportion declining
pal_fun <- RColorBrewer::brewer.pal(11, "PiYG") |>
  rev() |>
  colorRampPalette()
pal_decline <- pal_fun(length(decline_bins) - 1)
cols <- cut(trends_biomes_sf$decline_prop, decline_bins)
trends_biomes_sf$decline_prop_cols <- pal_decline[match(cols, levels(cols))]

for (b in unique(trends_biomes_sf$breeding_biome_label)) {
  pts <- filter(trends_biomes_sf, breeding_biome_label == b)

  # circles
  spp_bins <- quantile(pts$n_species, seq(0, 1, by = 0.05)) |>
    unique()
  midpoint_radius <- sqrt(spp_bins[-length(spp_bins)] + diff(spp_bins) / 2)
  spp_bins[1] <- 0
  spp_bins[length(spp_bins)] <- Inf
  radius_range <- c(2000, max_radius)
  circle_rad <- scales::rescale(midpoint_radius, to = radius_range) |>
    unname()
  rads <- cut(pts$n_species, spp_bins)
  pts$circle_radius <- circle_rad[match(rads, levels(rads))]
  circles <- st_buffer(pts, dist = pts$circle_radius)
  circles_sp <- st_buffer(pts, dist = max_radius)

  # make dot maps: trend
  path(figures_dir, glue("figure-04_{b}_abd-ppy.png")) |>
    png(width = 2400, height = 2400)
  par(mar = c(0, 0, 0, 0), bg = bg_col)

  # basemap
  plot(bb, col = NA, border = NA, axes = FALSE, bty = "n", reset = FALSE)
  plot(basemap$land,
       col = land_col, border = land_border_col, lwd = line_width,
       axes = FALSE, bty = "n", reset = FALSE, add = TRUE)
  plot(basemap$lakes,
       col = bg_col, border = land_border_col,
       lwd = line_width,
       axes = FALSE, bty = "n", reset = FALSE, add = TRUE)

  # trends
  plot(st_geometry(circles),
       col = circles$abd_ppy_cols, border = NA,
       axes = FALSE, bty = "n", reset = FALSE, add = TRUE)

  # borders
  plot(basemap$state_lines,
       col = line_col, lwd = line_width,
       axes = FALSE, bty = "n", reset = FALSE, add = TRUE)
  plot(basemap$country_lines,
       col = line_col, lwd = line_width * 2,
       axes = FALSE, bty = "n", reset = FALSE, add = TRUE)
  box()

  dev.off()
}

distinct(trends, breeding_biome, species_code, season) |>
  inner_join(ebird_taxonomy, by = "species_code") |>
  select(breeding_biome, species_code, common_name, season) |>
  arrange(breeding_biome, species_code) |>
  write_csv(path(figures_dir, "figure-04_trends-species_breeding-biomes.csv"))

# abd ppy legend
lbl_brks <- seq(-1, 1, length.out = length(ppy_bins))
lbls <- c(paste("<=", -4), 0, paste(">=", 4))
png(path(figures_dir, "figure-04_legend_abd-ppy.png"), width = 400, height = 1200)
par(mar = c(0, 0, 0, 0), bg = "#ffffff")
plot(bb, border = FALSE)
t <- "Abundance Trend 2007-2021 [% change / year]"
fields::image.plot(lbl_brks, breaks = lbl_brks, col = pal,
                   axes = FALSE, box = FALSE, legend.only = TRUE,
                   smallplot = c(0.5, 0.6, 0.05, 0.95),
                   axis.args = list(at = c(-1, 0, 1),
                                    labels = lbls,
                                    cex.axis = 2.8, lwd.ticks = 0),
                   legend.args = list(text = t, side = 2, cex = 3.5, line = 2))
dev.off()
