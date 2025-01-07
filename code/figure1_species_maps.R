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

data_dir <- "data"
outputs_dir <- "outputs"
figures_dir <- "figures"

# read in trends estimates with breeding biomes and srd information
trends <- path(data_dir, "ebird-trends_2021_srd-biomes.parquet") |>
  read_parquet()

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

# figure 1a: species maps at different scales ----

# srd cells within each bcr
bcr_lookup <- path(data_dir, "summary-regions_srd-lookup.parquet") |>
  read_parquet() |>
  filter(str_starts(region_code, "BCR"))
bcr_srd <- bcr_lookup |>
  group_by(srd_id) |>
  top_n(n = 1, wt = coverage_fraction) |>
  ungroup() |>
  distinct(region_code, srd_id) |>
  inner_join(distinct(trends, species_code, season,
                      srd_id, latitude, longitude),
             by = "srd_id", relationship = "many-to-many") |>
  select(species_code, season, region_code, srd_id, longitude, latitude)
species_srd <- bcr_srd |>
  distinct(species_code, season, srd_id, latitude, longitude)
problems <- bcr_srd |>
  group_by(species_code, srd_id) |>
  filter(n() > 1)
stopifnot(nrow(problems) == 0)

# trends aggregated to bcr
trends_bcr <- trends |>
  inner_join(bcr_lookup, by = "srd_id", relationship = "many-to-many") |>
  group_by(species_code, season, region_code) |>
  summarize(abd_ppy = sum(abd_ppy_median * abd * coverage_fraction) / sum(abd * coverage_fraction),
            .groups = "drop") |>
  inner_join(bcr_srd, by = c("species_code", "season", "region_code"),
             relationship = "many-to-many") |>
  mutate(scale = "bcr") |>
  select(scale, species_code, season, srd_id, longitude, latitude, abd_ppy)
problems <- trends_bcr |>
  group_by(species_code, srd_id) |>
  filter(n() > 1)
stopifnot(nrow(problems) == 0)

# rangewide trends
species_srd <- trends |>
  semi_join(trends_bcr, by = c("species_code", "season", "srd_id")) |>
  distinct(species_code, season, srd_id, longitude, latitude)
trends_rw <- trends |>
  semi_join(trends_bcr, by = c("species_code", "season", "srd_id")) |>
  group_by(species_code, season) |>
  summarize(abd_ppy = sum(abd * abd_ppy_median) / sum(abd), .groups = "drop") |>
  inner_join(species_srd, by = c("species_code", "season"),
             relationship = "one-to-many") |>
  mutate(scale = "rangewide") |>
  select(scale, species_code, season, srd_id, longitude, latitude, abd_ppy)
problems <- trends_rw |>
  group_by(species_code, srd_id) |>
  filter(n() > 1)
stopifnot(nrow(problems) == 0)

# combine
focal_species <- c("houfin", "grbher3", "wooduc", "houwre")
trends_scales <- trends |>
  semi_join(trends_bcr, by = c("species_code", "season", "srd_id")) |>
  mutate(scale = "srd") |>
  select(scale, species_code, season, srd_id, longitude, latitude,
         abd_ppy = abd_ppy_median) |>
  bind_rows(trends_bcr, trends_rw) |>
  filter(species_code %in% focal_species)
rm(trends_bcr, trends_rw)

# circles
trends_scales_circles <- trends_scales |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs = crs) |>
  st_buffer(dist = max_radius)

# ppy bins
ppy_bins <- seq(0, 4, by = 0.5)
ppy_bins[1] <- min(abs(trends_scales_circles$abd_ppy)) / 10
ppy_bins[length(ppy_bins)] <- max(abs(trends_scales_circles$abd_ppy)) * 10
ppy_bins <- c(-rev(sort(ppy_bins)), sort(ppy_bins))
pal <- ebirdst_palettes(length(ppy_bins) - 1, type = "trends")
cols <- cut(trends_scales_circles$abd_ppy, ppy_bins)
trends_scales_circles$abd_ppy_cols <- pal[match(cols, levels(cols))]

fig_dir <- path(figures_dir, "figure1")
dir_create(fig_dir)
for (species in focal_species) {
  for (s in c("srd", "bcr", "rangewide")) {
    ppy_circles <- filter(trends_scales_circles,
                          species_code == species,
                          scale == s)

    # make dot maps: ppy
    path(fig_dir, glue("{species}_abd-ppy_{s}.png")) |>
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
    plot(st_geometry(ppy_circles),
         col = ppy_circles$abd_ppy_cols, border = NA,
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
}
