This repository contains code and data for _North American bird declines are greatest where species are most abundant_ (Johnston et al. 2025).

## Abstract

Efforts to address declines of North American birds have been constrained by limited availability of fine-scale information about population change. Using participatory science data from eBird, we estimated continental population change and relative abundance at 27 km resolution for 495 bird species from 2007-2021. Results reveal high and previously undetected spatial heterogeneity in trends; although 75% of species were declining, 97% of species showed separate areas of significantly increasing and decreasing populations. Populations tended to decline most steeply in strongholds where species were most abundant, yet they fared better where species were least abundant. These high-resolution trends improve our ability to understand population dynamics, prioritize recovery efforts, and guide conservation at a time when action is urgently needed.

## Setup

To install all necessary packages to run the R scripts in the `code/` directory, run the following

```         
install.packages(readLines("requirements.txt"))
```

All R scripts in `code/` assume that the working directory is the top-level project directory, i.e. the one containing this README file. If you are using RStudio, open the `.Rproj` file and your working directory will be set correctly. If you're not using RStudio, you can set your working directory manually, e.g. using

```         
setwd("/Users/Alison/Documents/REPOS/2024-trends-abundance-paper/")
```

## Code

The R scripts included in the `code/` directory of this repository reproduce the figures and results from the paper.

- `figure-01_species-maps.R`: generates trend maps for Great Blue Heron, Wood Duck, and House Wren aggregated to three different spatial scales: range-wide, regional (Bird Conservation Region), and landscape scales (27 km × 27 km grid cells).
- `figure-02_fishbone.R`: generates a bar plot showing the minimum, maximum, median, and interquartile range of 27 km trend estimates for each of the 495 species in the study.
- `bam_linear_analysis.R`: Ali to provide a brief description of what this script does. This script generates a set of CSVs that are required by several of the subsequent scripts, therefore, it should be run prior to running the scripts used to generate Figures 3 and S14-S16.
- `figure-03_bam-abd.R`: Ali to provide a brief description of what this script does. Depends on `bam_linear_analysis.R`.
- `figure-04_community-maps.R`: generates community-level trend maps showing the mean trend across species breeding in six different biomes.
- `figure-S01_variation-within-biome.R`: generates a histogram showing the distribution of the range of 27 km trend estimates within each BCR-species combination.
- `figure-S02_trend-vs-prop-range-dec.R`: calculates the range-wide trend and proportion of 27 km grid cells that show a declining trend, then generates a boxplot showing the distribution of the proportion of declining cells binned by range-wide trend. TODO: missing the red box.
- `figure-S03-S04_density-prop-range-dec.R`: generates density plots showing the distribution of range-wide trend and the distribution of the proportion of cells with declining trends. Density plots are grouped by biome and (in Figure S4) split by breeding and non-breeding season trends.
- `figure-S09-S10_power-error.R`: generate plots showing the power and error versus trend magnitude (S9) and relative abundance (S10) for a set of 4 example species  Hermit Thrush, Hooded Warbler, Redhead, and Snow Goose.
- `figure-S11_power-error_all-species.R`: generates plots showing the power and error versus both the trend magnitude and relative abundance for all species.
- `figure-S12_power-error-histogram.R`: generates a histogram showing, for different levels of trend magnitude, what proportion of non-zero 27 km trends have directions that are correctly and incorrectly classified, respectively.
- `figure-S13_rangewide-error.R`: generates a histogram showing the distribution of range-wide error across species.
- `figure-S14-S15_range-edge.R`: generate plots showing the relationship between the trend estimates in 27 km cells and the distance of those cells to the species' range edge: a histogram of the slopes of the relationship for each species and the linear relationships by species grouped by biome. Depends on `bam_linear_analysis.R`.
- `figure-S16_compare-edge-abd.R`: based on the results of the BAM analysis, generates distributions of species-specific effect sizes for (a) distance to range edge and (b) relative abundance when both variables were included in the same model. Depends on `bam_linear_analysis.R`.
- `numbers-for-text.R`: produces a variety of numbers that appear in the text of the paper.

## Data

The datasets required to run the above scripts are in the `data/` directory. They are primarily tabular data in either CSV or Parquet format (Parquet files can be read into R using the `arrow` R package), and there is a single GeoPackage containing spatial data. The datasets and associated columns are as follows

**`master_species_list_495.csv`:** a list of all 495 species included in the study, with one row per species
- `species_code`: 6-letter eBird species code.
- `common_name`: English common name.
- `family`: Latin name of family.
- `breeding_biome`: breeding biome of the species.
- `season`: season for which the trend was estimated.
- `start_date/end_date`: dates defining the season boundaries.

**`ebird-trends_2007-2021.parquet:** trends estimates for all species 27 km resolution for the 2007-2021 time period. In addition to the trends, a set of species-specific simulation-based metrics capturing the ability of the model to correctly classify trends is included. These metrics are estimated across a suite of simulation scenarios for each grid cell grouped into bins according to trend magnitude.
- `species_code`: 6-letter eBird species code.
- `season`: season for which the trend was estimated.
- `breeding_biome/breeding_biome_label`: breeding biome of species.
- `srd_id`: unique integer ID for each 27 km grid cell.
- `longitude/latitude`: coordinates of 27 km grid cell centers.
- `abd`: relative abundance estimates for 2014, the midpoint of the 2007-2021 time period.
- `abd_ppy_median`: median percent per year trend in relative abundance. Note these are percentages expressed out of 100, e.g. a value of 5 expresses a 5% change in relative abundance per year.
- `abd_ppy_{lower/lower}`: the 10th and 90th percentile of percent per year trend in relative abundance.
- `abd_ppy_{sd/var}`: the variance and standard deviation of trend estimates based on a ensemble of 100 folds.
- `distance_to_edge_km`: distance between the grid cell and the range edge.
- `directional_power`: simulation based estimate of the proportion of trends estimates that are non-zero and correctly classify the trend direction.
- `directional_error_nz`: simulation based estimate of the proportion of non-zero trends estimates that incorrectly classify the trend direction.
- `bin_{min/max}`: trend magnitude bin boundaries over which the error and power metrics are estimated.

**`ebird-trends_abd-binned-performance_2021.parquet:** simulation-based power and error estimates for each species, grouped by relative abundances binned using ten equally spaced quantiles.
- `species_code`: 6-letter eBird species code.
- `season`: season for which the trend was estimated.
- `abd_bin_midpoint`: mid-point of the relative abundance bin.
- `n_nonzero`: number of non-zero trend estimates.
- `prop_nonzero`: proportion of trend estimates that are non-zero.
- `directional_power`: proportion of trends that are both non-zero and correctly classify the trend direction.
- `directional_error`: proportion of trends that are both non-zero and incorrectly classify the trend direction.
- `directional_power_nz`: proportion of non-zero trends that correctly classify the trend direction.
- `directional_error_nz`: proportion of non-zero trends that incorrectly classify the trend direction.

**`ebird-trends_ppy-binned-performance_2021.parquet:** simulation-based power and error estimates for each species, grouped into half percent trend magnitude bins.
- `species_code`: 6-letter eBird species code.
- `season`: season for which the trend was estimated.
- `ppy_bin_midpoint`: mid-point of the trend magnitude bins.
- `n_nonzero`: number of non-zero trend estimates.
- `prop_nonzero`: proportion of trend estimates that are non-zero.
- `directional_power`: proportion of trends that are both non-zero and correctly classify the trend direction.
- `directional_error`: proportion of trends that are both non-zero and incorrectly classify the trend direction.
- `directional_power_nz`: proportion of non-zero trends that correctly classify the trend direction.
- `directional_error_nz`: proportion of non-zero trends that incorrectly classify the trend direction.

**`ebird-trends_range-wide_folds_2021.parquet:** range-wide trend estimates for each of the 100 folds that make up the ensemble.
- `species_code`: 6-letter eBird species code.
- `season`: season for which the trend was estimated.
- `breeding_biome`: breeding biome of species.
- `fold`: integer indicating the fold number.
- `n_cells`: number of 27 km cells within the range.
- `abd_ppy`: abundance-weighted mean rangewide percent per year trend.
- `prop_decline_range`: proportion of 27 km cells with declining trends within range.

**`ebird-trends_simulations_focal-species_2021.parquet:** simulation results for a subset of four species: Hermit Thrush, Hooded Warbler, Redhead, and Snow Goose.
- `species_code`: 6-letter eBird species code.
- `season`: season for which the trend was estimated.
- `scenario_id`: unique integer ID for each simulation scenario (1-10).
- `realization`: each simulation scenario is repeated for 10 realization, this integer identifies the realization.
- `srd_id`: unique integer ID for each 27 km grid cell.
- `simulated`: simulated trend.
- `estimated`: median model estimated trend.
- `lower`: 10th percentile of model estimated trends.
- `upper`: 90th percentile of model estimated trends.
- `nonzero`: whether the confidence intervals overlap zero.

**`bcr-srd-lookup.parquet:** lookup table identifying which Bird Conservation Region (BCR) each 27 km grid cell falls into.
- `region_code`: unique ID for each BCR.
- `srd_id`: unique integer ID for each 27 km grid cell.
- `coverage_fraction`: proprtion of the 27 km grid cell covered by the BCR.

**`basemap.gpkg:** spatial data for creating basemaps in the mapping scripts. All data come from the [Natural Earth](https://www.naturalearthdata.com/) project. Layers included are:
- `land`: land boundary polygon.
- `country_lines`: country boundary lines.
- `state_lines`: state boundary lines.s
- `lakes`: lake polygons.
- `global_bb`: global boundary polygon.