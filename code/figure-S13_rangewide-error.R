library(arrow)
library(fs)
library(mgcv)
library(scales)
library(maps)
library(geosphere)

data_dir <- "data"
figures_dir <- path("figures", "figure-S13")
dir_create(figures_dir)

# power & error calculated as trend magnitude
D <- path(data_dir, "ebird-trends_2007-2021.parquet") |>
  read_parquet()
n_species <- length(unique(D$species_code))

# Order by range size
range.size <- table(D$species_code)
range.sizes <- order(range.size, decreasing =T)
species_code_by_rangesize <- names(range.size[range.sizes])

# Save linear abd effect summary
spp_trend_disp <- matrix(NA, n_species,2)
spp_trend_disp <- as.data.frame(spp_trend_disp)
names(spp_trend_disp) <-  c("sd", "iqr" )

expected_performance_bybin <- NULL
expected_performance <- NULL

for (iii.spp in seq_along(range.sizes)) {
  species_code <- names(range.size[range.sizes][iii.spp])
  message(paste(iii.spp, species_code))
  
  gam_data <- D[D$species_code == species_code, ]
  gam_data <- as.data.frame(gam_data)
  
  spp_expected_performance <- NULL
  for (iii.bin in seq_along(unique(gam_data$bin_min))){
    magnitude_bin_index <- 
      0.5*(iii.bin-1) <= abs(gam_data$abd_ppy_median) &  
      abs(gam_data$abd_ppy_median) < iii.bin*0.5
    if (sum(magnitude_bin_index)>0) {
      ttt <- gam_data[ magnitude_bin_index, ]
      ttt$directional_error_nz[ is.na(ttt$directional_error_nz) ] <- 0
      is_nz <- ttt$abd_ppy_lower < 0 & ttt$abd_ppy_upper > 0
      num_nz <- sum(magnitude_bin_index) - sum(is_nz)
      ttt.df <- data.frame(
        species_code = ttt$species_code[1],
        bin_min = ttt$bin_min[1],
        bin_max = ttt$bin_max[1],
        num_nz = num_nz, 
        dir_error_num = ttt$directional_error_nz[1] * num_nz, 
        total = sum(magnitude_bin_index) )
      
      spp_expected_performance <- rbind(spp_expected_performance, ttt.df)
    }
  } # iii.bin
  
  expected_performance_bybin <- rbind(
    expected_performance_bybin, spp_expected_performance)
  ttt.df <- data.frame(
    species_code = spp_expected_performance$species_code[1],
    dir_error = sum(spp_expected_performance$dir_error_num) / 
      sum(spp_expected_performance$total) )
  expected_performance <- rbind(expected_performance, ttt.df) 
  
} # iii.spp


plot_loc <- path(figures_dir, paste0("figure-S13_rangewide-error.tif"))
tiff(plot_loc, width = 8, height = 8, units = "cm", pointsize = 9, res = 600)
nbins <- 25
hist(expected_performance$dir_error*100,  
     breaks = nbins,
     xlab = "Expected Range-wide Error(%)",
     main = "") 
abline(v=0)
dev.off()

mean(expected_performance$dir_error*100)
sd(expected_performance$dir_error*100)
summary(expected_performance$dir_error*100)
