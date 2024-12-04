#########################################################
## fit the linear mixed model for abundance and 
## create figure 3

library(tidyverse)
library(arrow)
library(mgcv)
library(scales)
library(maps)
library(geosphere)
library(fs)

base_dir <- "/Users/Alison/Documents/REPOS/2024_trends-abundance-paper/"
data_dir <- path(base_dir, "data")
outputs_dir <- path(base_dir, "outputs")
figures_dir <- path(base_dir, "figures/figure3/bam")
dir.create(figures_dir, recursive = TRUE)

#########################################################
## load data

# # species lookup
species <- read_csv(path(data_dir, "/master_species_list_495.csv"), na = "") |>
  select(species_code, common_name, breeding_biome)

trends <- path(data_dir, "ebird-trends_2021_weights.parquet") |>
  read_parquet() |>
  mutate(log10_abd = log10(abd),
         weight = 1 / abd_ppy_var) |>
  mutate(log10_distance_to_edge_km = log10(distance_to_edge_km)) |>
         # breeding_biome = factor(breeding_biome),
         # species_code = factor(species_code),
  left_join(species, by = "species_code")

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
## FUNCTIONS TO RUN THE BAM MODELS FOR TREND AND ABUNDANCE



# -------------------------------------------------------
# create a number that is the number of knots for the 
# spatial component of the models
# this is a function of the range size of the species
# with a min and max cap

create_k_spatial <- function(mod_data, gam_min_ss = 100, gam_max_ss = 10000){
	# set spatial knots as function of range size
	k.spatial <- round((nrow(mod_data) - gam_min_ss)/gam_max_ss * 175) + 25

	return(k.spatial)
}




# -------------------------------------------------------
# a function that fits one single model. it takes as 
# parameter inputs the data and the 'model_type'
# model_type is one of: lm_linear, bam_linear, gamlss_linear, gamlss_smooth
# the first part of the model_type is the model function, the second
# part is the type of fit for the log10_abd parameter. 

fit_1_mod <- function(mod_data, rho, 
	k.spatial = NULL, gam_min_ss = 100, gam_max_ss = 10000, 
	k.base = 3, 
	species_code = NULL, 
	model_type = "bam_linear", 
	covariate = "log10_abd",
	...){

	if(is.null(species_code) && any(colnames(mod_data)=="species_code")) species_code <- mod_data$species_code[1]

	# reduce data to required columns
	#	mod_data <- mod_data |> dplyr::select(abd_ppy_median, abd_ppy_sd, log10_abd, directional_error_nz, longitude, latitude, species_code)

	# add covariate column
	if(covariate == "log10_abd") mod_data$covariate <- mod_data$log10_abd
	if(covariate == "log10_distance_to_edge_km") mod_data$covariate <- mod_data$log10_distance_to_edge_km

	# set spatial knots as function of range size
	if(is.null(k.spatial)){
			k.spatial <- create_k_spatial(mod_data = mod_data, gam_min_ss = gam_min_ss,  gam_max_ss = gam_max_ss)
	}

	# fit bam model
	if(model_type == "bam_linear"){

		b_mod <- bam(abd_ppy_median ~ covariate +
							s(longitude, latitude, k=k.spatial, 
								bs="gp", m=c(1, rho)),
					weights = abd_ppy_sd / mean(abd_ppy_sd),
					data = mod_data, 
					discrete = TRUE)
	}

	if(model_type == "gamlss_linear"){
		b_mod <- gam(list(
			# Mean
				abd_ppy_median ~  
					covariate +
					s(longitude, latitude, k=k.spatial, 
						bs="gp", 
						m=c(1, rho)),
			# Variance
				~ s(abd_ppy_sd, k=k.base, bs="ds", m=c(1,0)) + 
					s(directional_error_nz, k=k.base, bs="ds", m=c(1,0)) ), 

			data = mod_data,
			family = gaulss())
	}

	if(model_type == "gamlss_smooth"){
		b_mod <- gam(list(
			# Mean
				abd_ppy_median ~  
					s(covariate,
						k=k.base, 
						bs="ds", m=c(1,0)) +
					s(longitude, latitude, k=k.spatial, 
						bs="gp", 
						m=c(1, rho)),
			# Variance
				~ s(abd_ppy_sd, k=k.base, bs="ds", m=c(1,0)) + 
					s(directional_error_nz, k=k.base, bs="ds", m=c(1,0) )), 
			data = mod_data,
			family = gaulss())

	}

	if(model_type == "lm_linear"){
		b_mod <- lm(abd_ppy_median ~ covariate,
					weights = abd_ppy_sd / mean(abd_ppy_sd),
					data = mod_data)
	}

	# return required metric
	return(b_mod)
}



# -------------------------------------------------------
# a function to create a sequence of rho values to test 
# different models with different spatial correlation sturctures
# it depends on the range size of the species and the number of 
# rho parameters requested. 

create_rho_sequence <- function(mod_data, length_rho_seq = 10){

	# reduce data to required columns
#	mod_data <- mod_data |> dplyr::select(abd_ppy_median, abd_ppy_sd, log10_abd, directional_error_nz, longitude, latitude)

	# Max SS (for computational complexity)
	rho_cv_max_ss <- 2000

	# Construct spatial correlation range parameter grid
	# based on empricical quantiles of interpoint distances
	# approx distance bet points in A and B in meters
	A <- data.frame(
		  lat = mod_data$latitude,
		  lon = mod_data$longitude )
	A <- B <- A[ complete.cases(A), ] |> distinct()

	dist_matrix <- distm(
			cbind(B$lon, B$lat), 
			cbind(A$lon, A$lat), 
			fun = distHaversine)

	# Convert to km & remove zero distances
	dist_matrix_km <- dist_matrix/1000 
	dist_matrix_km[ dist_matrix_km ==0 ] <- NA

	k <- 3
	rho_sequence <- quantile(dist_matrix_km, 
		probs = seq(
				from = (0.01)^(1/k), 
				to = (0.9)^(1/k), 
				length=length_rho_seq)^k, 
			na.rm=TRUE)
	# Set first distance at 30km so smallest range
	# includes nearest neighbors		

	# Set first distance at 30km so smallest range
	# includes nearest neighbors, but is not less than 30km
	big_index <- rho_sequence > 30
	rho_sequence <- rho_sequence[big_index]
	if (rho_sequence[1] > 30) rho_sequence <- c(30, rho_sequence)

	# Translate back to a very approximate distance based on 
	# lat & lon as cartesian coordinates. 
	# The simplest approximation is 1 deg lon change is ~ 111km
	rho_sequence_ll <- rho_sequence / 111

	return(rho_sequence_ll)
}



# -------------------------------------------------------
# a function that finds the optimal rho for a dataset
# and model_type. outputs two 'optimal' rhos and 
# the dataset of REML and rho values tested. 

find_opt_rho <- function(mod_data, length_rho_seq = 10, 
	k.spatial = NULL, gam_min_ss = 100, gam_max_ss = 10000, 
	species_code = NULL, 
	plots_dir, plot_flag = FALSE, 
	model_type = "bam_linear",
	rho_method = "first_within_2sd",
	covariate = "log10_abd",
	...){

	if(is.null(species_code)) species_code <- mod_data$species_code[1]

	# reduce data to required columns
#	mod_data <- mod_data |> dplyr::select(abd_ppy_median, abd_ppy_sd, log10_abd, directional_error_nz, longitude, latitude)

	# create rho sequence
	rho_sequence_ll <- create_rho_sequence(mod_data, length_rho_seq = length_rho_seq)

	if(is.null(k.spatial)){
		k.spatial <- create_k_spatial(mod_data, gam_min_ss, gam_max_ss)
	}

	# run the bam for each value of rho
	REML <- rep(NA, length(rho_sequence_ll))

	for(r in 1:length(rho_sequence_ll)) {
		rho <- rho_sequence_ll[r]
		mod <- fit_1_mod(mod_data, rho = rho, k.spatial = k.spatial, 
			model_type = model_type, covariate = covariate)
		REML[r] <- mod$gcv.ubre
	}

	# combine rho and REML into data frame
	df <- data.frame(rho_km = rho_sequence_ll*111, REML)

	# absolute minimum REML
	min_rho_km <- df$rho_km[ which.min(df$REML) ]
	min_reml <- min(df$REML)

	# -------------------------------------------------------
	# option 1: rho_km with the minimum REML
	if(rho_method %in% c("first_min")){
		final_rho <- min_rho_km
	}

	# -------------------------------------------------------
	# option 4
	# first point at which it is within 2*sd of the mean of the last third of the time series
	# (nb: needs min 20 for this to be reasonable)
	
	if(rho_method %in% c("first_within_2sd")){

		# find mean and sd of stable part of time series
		mn_stable <- mean(df$REML[floor(nrow(df) - nrow(df)/3):nrow(df)])
		sd_stable <- sd(df$REML[floor(nrow(df) - nrow(df)/3):nrow(df)])

		w4 <- df$REML < (mn_stable + sd_stable*2)
		first_within_2sd <- min(which(w4))
		final_rho <- df$rho_km[first_within_2sd]
	}

	# Save Diagnostic plot and the rho estimates
	if(plot_flag){
		png(filename = 
			path(plots_dir,
				paste0(model_type, "_", species_code, "_rho.cv.diag.png")),
			width = 700, height = 700 )

			par(mfrow=c(1,1), cex=1.25, mar=c(6,6,6,6))		
		
			plot(df$rho_km, df$REML, 
				xlab = "rho (km)",
				ylab = "REML", 
				type = "o", 
				pch = 16,
				main = species_code)

			abline(v = min_rho_km, lwd=3, col="blue")
			abline(v = final_rho, lwd=2, col="red")

		dev.off()
	}

	ret <- list(final_rho = final_rho, df = df)
	return(ret)
}



# -------------------------------------------------------
# a wrapper function that first finds an optimal rho for 
# a dataset and model_type. then runs and outputs a 
# final model using the optimal rho. 
# currently hard-coded to use 'first_low_rho_km' 
# 'min_rho_km' is an alternative option. 
# (see function find_opt_rho for details)

fit_opt_model <- function(mod_data, species_code = NULL, 
				length_rho_seq = length_rho_seq, gam_min_ss = 100, gam_max_ss = 10000, 
				plots_dir, plot_rho_diagnostic = FALSE,
				model_type = "bam_linear", 
				rho_method = "first_within_2sd",
				covariate = "log10_abd",
				...){

	if(is.null(species_code)) species_code <- mod_data$species_code[1]

	# subset to required columns
	#	mod_data <- mod_data |> dplyr::select(abd_ppy_median, abd_ppy_sd, log10_abd, directional_error_nz, longitude, latitude, species_code)

	# find spatial knot number
	k.spatial <- create_k_spatial(mod_data, gam_min_ss = gam_min_ss, gam_max_ss = gam_max_ss)

	# find optimal rho (if using a spatial model type)
	# (rho irrelevant if linear model)
	if(model_type != "lm_linear"){
		opt_rho <- find_opt_rho(mod_data, length_rho_seq = length_rho_seq, k.spatial = k.spatial, 
			species_code = species_code, 
			plots_dir = plots_dir, plot_flag = plot_rho_diagnostic, 
			model_type = model_type, 
		    rho_method = rho_method,
		    covariate = covariate,
			...)
		rho_use <- opt_rho$final_rho
	}
	if(model_type == "lm_linear") rho_use <- NA 

	# run final model
	final_model <- fit_1_mod(mod_data, rho = rho_use, k.spatial = k.spatial, model_type = model_type,
		species_code = species_code, covariate = covariate)

	ret <- list(final_model = final_model, model_type = model_type)
	return(ret)

}




# -------------------------------------------------------
# function to extract the log10_abd coefficient from the model if linear model
# and the edf if a smooth. also keep the test statistic (T/Chisq) 
# and p-value

extract_effect <- function(model_output, variable = "covariate"){

	variable_name <- variable
	if(variable == "intercept") variable_name <- "(Intercept)"

	model_type <- model_output[["model_type"]]
	model_obj <- model_output[["final_model"]]

	if(model_type == "lm_linear"){
		tab <- summary(model_obj)$coefficients
		res <- tab[rownames(tab) == variable_name,] |> as.numeric()
	}

	if(model_type %in%  c("gamlss_linear", "bam_linear")){
		tab <- summary(model_obj)$p.table
		res <- tab[rownames(tab) == variable_name,] |> as.numeric()
	}

	if(model_type == "gamlss_smooth"){
		if(variable == "covariate"){
			tab <- summary(model_obj)$s.table
			res <- tab[rownames(tab) == paste0("s(", variable_name, ")"),] |> as.numeric()
		} 
		if(variable == "intercept"){
			tab <- summary(model_obj)$p.table
			res <- tab[rownames(tab) == variable_name,] |> as.numeric()
		} 
	}

	names(res) <- colnames(tab)
	return(res)

}


# -------------------------------------------------------
# function to plot the species scatterplot of data, sized 
# by the weights and with the fitted effect of log10_abd

plot_model <- function(mod_data, model_output, species_code = NULL, species_name = NULL, plots_dir, covariate = "log10_abd"){

	if(is.null(species_code)) species_code <- mod_data$species_code[1]

	# collect species name for plots
	spec_code <- species_code
	if(is.null(species_name) & "common_name" %in% colnames(species)){
		species_name <- species |> filter(species_code == spec_code) |> pull(common_name)
		species_name <- gsub("\'", "", species_name, fixed = TRUE)
		species_name <- gsub("\ ", "_", species_name, fixed = TRUE)
	}

	# add covariate column
	if(covariate == "log10_abd") mod_data$covariate <- mod_data$log10_abd
	if(covariate == "log10_distance_to_edge_km") mod_data$covariate <- mod_data$log10_distance_to_edge_km

	# subset to required columns
	# mod_data <- mod_data |> dplyr::select(abd_ppy_median, abd_ppy_sd, log10_abd, directional_error_nz, longitude, latitude, species_code)

	model_type <- model_output[["model_type"]]
	model_obj <- model_output[["final_model"]]

	lo <- ifelse(model_type == "gamlss_smooth", 1000, 10)
	nd <- data.frame(covariate = seq(min(mod_data$covariate), max(mod_data$covariate), length.out = lo))

	# plot scatterplot of abundance vs trend
	# set up plot (to hopefully make pretty for many species)
	xlimits <- range(mod_data$covariate)
	la_vals <- seq(-5, 5, by = 1)

	# if(xlimits[2] - xlimits[1])
	# la_vals <- seq(-6, 6, by = 2)

	ylimits <- range(mod_data$abd_ppy_median)
	if(ylimits[2]< -10) ylimits[2] <- -10
	if(ylimits[1]> 10) ylimits[1] <- 10
	ylimits <- c(-10, 10)


	# extract coefficient
	coef <- extract_effect(model_output, variable = "covariate")[1]
	int <- extract_effect(model_output, variable = "intercept")[1]

	# create useful plot name
	formatted_slope <- paste0(ifelse(model_type == "gamlss_smooth", "s", ifelse(coef<0, "n", "p")), "_",
					format(round(abs(coef), digits = 3), nsmall = 3))
	species_plot_name <- paste0(formatted_slope, "_", species_code, "_", species_name)
	plot_loc <- path(plots_dir, paste0("species_scatter_", species_plot_name, ".png"))

	# create file image and set up plot
	png(plot_loc, width = 8, height = 8, units = "cm", pointsize = 9, res = 600)
	par(mar = c(5, 5, 3, 3))

		# create plot
		plot(0, 0, col = "white", 
		    xlim = xlimits, 
		    xaxt = "n", xlab = "Relative abundance",
		    ylim = ylimits,
		    yaxt = "n", ylab = "Population trend", 
		    main = species_name, 
		    mgp = c(1.5, 0.7, 0))
		axis(side = 1, at = la_vals, labels = 10^la_vals, tcl = -0.3, mgp = c(1.5, 0.5, 0))
		axis(side = 2, at = c(-10, 0, 10), las = 1, tcl = -0.3, mgp = c(1.5, 0.7, 0))

		abline(h = 0, col = alpha("firebrick", 0.2), lwd = 2)

		points(mod_data$covariate, mod_data$abd_ppy_median, 
	        pch = 16, col = alpha("steelblue", 0.3),
			cex = sqrt(mod_data$abd_ppy_sd / mean(mod_data$abd_ppy_sd)))

		abline(int, coef, col = "steelblue", lwd= 1.5)

	dev.off()

}


# -------------------------------------------------------
# wrapper funciton for everything together 

run_whole_thing_per_species <- function(mod_data, model_type = "lm_linear", 
	plots_dir_rho, plot_rho_diagnostic = TRUE,
	plots_dir_scatter, plot_scatter = TRUE,
	length_rho_seq = 10,
	file_suffix = "",
	gam_min_ss = 100, gam_max_ss = 10000,
	covariate = "log10_abd",
	...){

	time1 <- Sys.time()

	species_code <- mod_data$species_code[1]
	print(paste("running for species", species_code))
	print(paste("number pixels =", nrow(mod_data)))

	# collect species name for plots
	spec_code <- species_code
	species_name <- species_code

	# this code for pretty species names works on examples, but currently
	# buggy in bulk for reasons I can't figure out. 
	# if("common_name" %in% colnames(species)){
	# 	species_name <- species |> filter(species_code == spec_code) |> pull(common_name)
	# 	species_name <- gsub("\'", "", species_name, fixed = TRUE)
	# 	species_name <- gsub("\ ", "_", species_name, fixed = TRUE)
	# } 

	# check if species results are already on file
	results_file <- path(outputs_dir, paste0("species_coefs_", model_type, "_", file_suffix, ".csv"))
	species_results_exist <- FALSE
	if(file.exists(results_file)){
		check <- read_csv(results_file)
		species_results_exist <- ifelse(any(check$species_code == spec_code), TRUE, FALSE)	
	}

	if(species_results_exist) { print("species results already on file"); coef_summary <- NULL }

	if(species_results_exist == FALSE){

		# find optimal rho and fit model
		fit_mod <- fit_opt_model(mod_data = mod_data, model_type = model_type, 
			length_rho_seq = length_rho_seq, 
			plots_dir = plots_dir_rho, plot_rho_diagnostic = plot_rho_diagnostic,
			gam_min_ss = gam_min_ss, gam_max_ss = gam_max_ss, covariate = covariate, ...)

		# plot model
		if(plot_scatter){
			plot_model(mod_data = mod_data, model_output = fit_mod, 
				species_code = species_code, species_name = species_name, 
				plots_dir = plots_dir_scatter, covariate = covariate, ...)
		}

		# save slope estimates
		coefs <- c(extract_effect(fit_mod, variable = "covariate"))
		coef_df <- data.frame(est = coefs[1], se = coefs[2], t_val = coefs[3], p_val = coefs[4], species_code = species_code)

		# append to file
		if(file.exists(results_file)) write_csv(coef_df, results_file, append = TRUE)
		if(!file.exists(results_file)) write_csv(coef_df, results_file, col_names = TRUE)

		time2 <- Sys.time()

		print(paste("finished species", species_code))
		print(time2 - time1)

		return(coef_df)
	}
}




#########################################################
## RUN FOR ALL SPECIES 

sample_size <- data.frame(species_code = names(table(trends$species_code)), n_pixels = as.numeric(table(trends$species_code)))
trends_n <- trends |> 
		left_join(sample_size)

# subset to do quick test before full run if useful 
trends_sub <- trends_n |>
		filter(species_code %in% c("abetow", "acafly", "acowoo", "aldfly", "allhum", "altori"))


run_name <- "log10_abd_40"
run_name <- "testtest"

scatterplots_dir <- path(outputs_dir,
 			"scatterplots", run_name)
if(!file.exists(scatterplots_dir)) dir.create(scatterplots_dir)

opt_rho_plots_dir <- path(outputs_dir,
 			"optimise_rho_plots", run_name)
if(!file.exists(opt_rho_plots_dir)) dir.create(opt_rho_plots_dir)

results_bam_linear <- trends_sub |>
		arrange(n_pixels) |>
		filter(n_pixels > 100) |>
		group_by(n_pixels, species_code) |>
		group_split() |>
		map_dfr(run_whole_thing_per_species, model_type = "bam_linear", 
		plots_dir_rho = opt_rho_plots_dir, plots_dir_scatter = scatterplots_dir, 
		length_rho_seq = 40, file_suffix = run_name, gam_min_ss = 100, 
		covariate = "log10_abd")


# change covariate = "log10_abd" to covariate = "log10_distance_to_edge_km"
# and rerun. 








# mod_data <- trends |> filter(species_code == "mexjay4")
# test <- run_whole_thing_per_species(mod_data, model_type = "bam_linear", 
# 		plots_dir_rho = opt_rho_plots_dir, plots_dir_scatter = scatterplots_dir, 
# 		length_rho_seq = 20, file_suffix = run_name, gam_min_ss = 100, 
# 		covariate = "log10_distance_to_edge_km")


# # -------------------------------------------------------
# # test run for one species

# run_name <- "test_edge_v2"

# scatterplots_dir <- path(outputs_dir,
#  			"scatterplots", run_name)
# if(!file.exists(scatterplots_dir)) dir.create(scatterplots_dir)

# opt_rho_plots_dir <- path(outputs_dir,
#  			"optimise_rho_plots", run_name)
# if(!file.exists(opt_rho_plots_dir)) dir.create(opt_rho_plots_dir)


# # run_whole_thing_per_species
# mod_data <- trends |> filter(species_code == "mexjay4")
# model_type <- "bam_linear"
# plots_dir_rho <- opt_rho_plots_dir
# plots_dir_scatter <- scatterplots_dir
# length_rho_seq <- 30
# file_suffix <- run_name
# covariate = "log10_distance_to_edge_km"

# #
# plots_dir <- plots_dir_rho
# plot_rho_diagnostic <- TRUE
# gam_min_ss <- 100
# gam_max_ss <- 10000

# # fit_opt_model
# species_code <- mod_data$species_code[1]
# gam_min_ss = 100; gam_max_ss = 10000
# rho_method = "first_within_2sd"

# # find_opt_rho
# plot_flag = TRUE


