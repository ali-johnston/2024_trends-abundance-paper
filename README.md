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

