# Author: Sergio M. Vicente-Serrano  <http://www.ipe.csic.es/vicente-serrano-s.m.>; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC <http://www.ipe.csic.es/hidrologia-ambiental>. Fergus Reig Gracia <http://fergusreig.es/>; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC <http://www.ipe.csic.es/hidrologia-ambiental/>
# Version: 1.0

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/> <http://www.gnu.org/licenses/gpl.txt/>.
#####################################################################


#' Performs quality control
#' Quality control: Stations with less than 20 years of data are removed, and using the 10 most correlated stations within 200 km, data with a percentile difference greater than 0.6 are discarded.
#'
#' @param data path to the data file
#' @param max_dist_eval maximum distance between two stations to use one to evaluate or complete the other
#'
#' @return data and coordinates with data that do not pass the control removed
#' @export
#'
mediterranean_calculations <- function(data, max_dist_eval){

	# Remove data if there are 5 or more consecutive months of zeros, if any of the months involved has less than 70% zeros
	delete_zero_data <- delete_zero(data = data$data)

	# Quality control, complete series and invalid loose data are removed
	control_data <- quality_control(data = delete_zero_data, coor = data$coor, max_dist = max_dist_eval, max_diff_anomaly = max_diff_anomaly, max_diff_anomaly_0 = max_diff_anomaly_0)

	return(control_data)
}

#' Performs a second fill
#' For each dataset (1870, 1900...)
#' Stations with more than 90 or 95% of filled data
#' Sort the stations by correlation (minimum 0.5) and fill using the 10 methods...
#' Stations without total fill are discarded
#'
#' @param file_data path to the data file
#' @param fillable_years years that can be filled with the station's monthly average
#' @param max_dist maximum allowed distance for filling
#'
#' @return None
#' @export
#'
second_data_fill_data <- function(file_data, fillable_years = 36, max_dist = NA){
	# Remove stations without 95% of data
	no_nas <- (dim(file_data$data)[1] - apply(file_data$data, c(2), sum_no_nas)) < (dim(file_data$data)[1] * (100 - percentage_filled_data) / 100)
	control_data <- list(data = file_data$data[, no_nas, drop = FALSE], coor = file_data$coor[no_nas, , drop = FALSE])

	# Second fill
	fill_data <- fill_series(control_data = control_data, min_correlation = min_second_correlation, max_dist = max_dist)

	# Remove incomplete stations
	no_nas <- dim(fill_data$data)[1] - apply(fill_data$data, c(2), sum_no_nas) == 0
	all_data <- list(data = fill_data$data[, no_nas, drop = FALSE], coor = fill_data$coor[no_nas, , drop = FALSE])

	# If there is no data, try to fill some series with itself
	if(dim(all_data$coor)[1] <= 6){
		# Fill 2 or 3 years of data with the station's own monthly average
		fill_data <- fill_unfillable_station(data = fill_data, fillable_years = fillable_years)

		# Remove incomplete stations
		no_nas <- dim(fill_data$data)[1] - apply(fill_data$data, c(2), sum_no_nas) == 0
		all_data <- list(data = fill_data$data[, no_nas, drop = FALSE], coor = fill_data$coor[no_nas, , drop = FALSE])
	}

	return(all_data)
}

#' Performs a second fill
#' For each dataset (1870, 1900...)
#' Stations with more than 90 or 95% of filled data
#' Stations are sorted by correlation (minimum 0.5), and filling is done using the 10 methods...
#' Stations without total fill are discarded
#'
#' @param data data and coordinates
#' @param max_dist_eval maximum distance for filling
#'
#' @return data and coordinates
#' @export
#'
second_data_fill <- function(data, max_dist_eval = NA){
  all_data <- list()
  i_ini <- names(data)[length(data)]
  for(i_ini in names(data)){
  	file_data <- data[[i_ini]]

  	# Years that can be filled with the station's monthly average
  	if(i_ini == paste0("start_", i_inis[length(i_inis)])){
  		fillable_years <- 24 # 2 years
  	}else{
  		fillable_years <- 36 # 3 years
  	}

  	data_fill <- second_data_fill_data(file_data = file_data, fillable_years = fillable_years, max_dist = max_dist_eval)
  	if(sum(!is.na(data_fill)) > 0){
  		all_data[[i_ini]] <- data_fill
  	}
	}
	return(all_data)
}

#' Homogenization - Alexanderson test
#' Reference series, compares and corrects
#' For each dataset (1870, 1900...)
#' The 5 most correlated series are chosen using the difference series
#' A weighted average is made with the 5 series (correlation * data1 + ...) / sum(correlations), which will be the reference series
#' Alexanderson will give a breakpoint and a ratio value to multiply the older part by... iterate while breakpoints are given
#' Save statistics on inhomogeneities. Basically, the number of data points changed in each series and the time of the inhomogeneity.
#' CSV with the number of data points changed and CSV with inhomogeneity point - all x 12 months
#'
#' @param file_data path to the data file
#' @param no_use_series series that will not be homogenized
#'
#' @return None
#' @export
#'
alexanderson_homogenize_data <- function(file_data, no_use_series = c()){

	max_i_snht <- 20
	significance_level <- 100 * (1 - 0.01)
	data_save <- file_data$data

	# Search for inhomogeneities only in series that were not already in previous periods
	all_series <- colnames(file_data$data)
	use_series <- all_series[!all_series %in% no_use_series]

	# Differences from one day to the next
	differences_0 <- file_data$data
	differences_1 <- rbind(rep(0, dim(differences_0)[2]), differences_0[c(1:(dim(differences_0)[1]-1)), ])
	differences <- differences_0 - differences_1

 	data_cor_list <- near_correlations(data = differences, coor = file_data$coor, max_dist = NA)
	months <- names(data_cor_list)

	break_points <- array(NA, dim = c(length(months), length(all_series), max_i_snht), dimnames = list(months, all_series, c(1:max_i_snht)))
	change_values <- array(NA, dim = c(length(months), length(all_series), max_i_snht), dimnames = list(months, all_series, c(1:max_i_snht)))

	i_month <- months[5]
	for(i_month in months){
	    data_cor <- data_cor_list[[i_month]]
    	date_month <- grepl(i_month, rownames(file_data$data))
    	data_save_month <- file_data$data[date_month, ]
    	data_month <- NA
    	i_snht <- 1
		while(i_snht < max_i_snht & (length(data_month) == 1 | sum(data_save_month != data_month, na.rm = TRUE) > 0)){
			# print(paste(i_month))
			data_month <- data_save_month
			i_series <- use_series[12]
			for(i_series in use_series){
	    		i_data_cor <- data_cor[, i_series]
	    		i_data_cor <- i_data_cor[!is.na(i_data_cor)]
	    		if(length(i_data_cor) > 0){
		    		i_estations <- min(length(i_data_cor), 5)
		    		i_order <- names(i_data_cor)[order(i_data_cor, decreasing = TRUE)][c(1:i_estations)] #Cinco estaciones más correlacionadas
		    		i_data_cor_select <- i_data_cor[i_order]

		    		# Weighted average
		    		average <- t(i_data_cor_select * t(data_month[, i_order])) / sum(i_data_cor_select)
		    		average <- apply(average, c(1), sum)

		    		# snht
		    		Q <- (data_month[, i_series] + 1) / (average + 1)
		    		Q[Q == Inf] = NA
		    		Q[Q == -Inf] = NA
		    		Q[is.nan(Q)] = NA		    		

		    		if(sum(!is.na(Q)) > 0 & length(unique(Q)) > 1){
			    		q_snht <- snht(cbind(1:length(Q), Q), significance_level) # https://sites.google.com/site/phaenggi23/rscripts
			    		break_serie <- q_snht$T0x
	      				pre_post_mean <- round(mean(data_month[c((break_serie + 1):(dim(data_month)[1])), i_series]) / mean(data_month[c(1:break_serie), i_series]), digits = 15)
	      				pre_post_mean[pre_post_mean == Inf] = NA
			    		pre_post_mean[pre_post_mean == -Inf] = NA
			    		pre_post_mean[is.nan(pre_post_mean)] = NA
			    		if (q_snht$T0 > q_snht$Tc && !is.na(pre_post_mean) && pre_post_mean != 1 && pre_post_mean != 0) {       			
		        			data_save_month[c(1:break_serie), i_series] <- pre_post_mean * data_month[c(1:break_serie), i_series]
		        			# print(paste(i_month, i_series, break_serie))
		        			break_points[i_month, i_series, i_snht] <- rownames(data_month)[break_serie]
		        			change_values[i_month, i_series, i_snht] <- break_serie
	    				} #if
					} #if
				} #if
			} #for
				i_snht <- i_snht + 1
		} #while
		data_save[date_month, ] <- data_save_month
	}#for
	return(list(data = data_save, break_points = break_points, change_values = change_values))
}

#' Alexanderson test for all available files (which have successfully passed the second fill)
#'
#' @param data data and coordinates
#' @param folder directory to save the output files
#'
#' @return data and coordinates
#' @export
#'
alexanderson_homogenize <- function(data, folder){
	dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  
  pre_stations <- NA
  file_data_pre <- NULL
  data_save <- list()
	i_ini <- names(data)[1]
	for(i_ini in names(data)){
    	file_data <- data[[i_ini]]
    	if(length(file_data) > 1 && dim(file_data$coor)[1] > 0){
    		if(dim(file_data$coor)[1] > 1){
	    		# The stations already in the previous group must be exactly the same
	    		if(!is.null(file_data_pre)){
	    			delete_stations <- pre_stations[!pre_stations %in% colnames(file_data$data)]
	    			if(length(delete_stations) > 0){
	    				file_data$data <- cbind(file_data$data, file_data_pre$data[rownames(file_data$data), delete_stations, drop = FALSE])
	    				file_data$coor <- rbind(file_data$coor, file_data_pre$coor[delete_stations, ])
	    			}
	    			file_data$data[, pre_stations] <- file_data_pre$data[rownames(file_data$data), pre_stations]
	    		}
	    		homogenize_data <- alexanderson_homogenize_data(file_data = file_data, no_use_series = pre_stations)
		    	
	    		# Generate pre-homogenization statistics
	    		mkTrend_pval <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_pval)

	  			mkTrend_slp <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_slp)

					percentage <- calc_data_year_month_station(data = file_data$data, calc_function = calc_percentage)

		    	data_save[[i_ini]] <- list(data = homogenize_data$data, coor = file_data$coor, change_values = homogenize_data$change_values, break_points = homogenize_data$break_points, mkTrend_pval_pre_homogenize = mkTrend_pval, mkTrend_slp_homogenize = mkTrend_slp, percentage_homogenize = percentage)
	    	}else{
	    		data_save[[i_ini]] <- list(data = file_data$data, coor = file_data$coor, change_values = NA, break_points = NA, mkTrend_pval_pre_homogenize = NA, mkTrend_slp_homogenize = NA, percentage_homogenize = NA)
	    	}
	    	pre_stations <- rownames(data_save[[i_ini]]$coor)
	    	file_data_pre <- data_save[[i_ini]]
				save_csvs(i_ini, folder_name = folder, data_save = data_save[[i_ini]]$data, coor_save = data_save[[i_ini]]$coor)
  		}#if
	}#for
	return(data_save)
}

#' Calculates statistics for the data
#' Monthly, seasonal, and annual trends, Trend package, sens.slope function
#' Add 1 to everything to avoid zeros
#' Significance, modifiedmk package, bbsmk function
#' National average series
#' SPI at scales 3, 12, and 24 for each series, important to ensure operations cannot be reversed
#'
#' @param file_data data and coordinates
#' @param data_ori original data
#'
#' @return None
#' @export
#'
calculate_statistics_data <- function(file_data, data_ori){

	mobile_trends_calc <- c(1871, 1901, 1931)
	i_ini <- read_years(rownames(file_data$data))[1]
	if(i_ini %in% mobile_trends_calc){ # Calculating mobile_trends_data is very slow, so it will only be done for some cases
		calc_mobile_trends_data = TRUE
	}else{
		calc_mobile_trends_data = FALSE
	}

	# National average series
	average <- apply(file_data$data, c(1), mean)

	# SPI at scales 3, 9, 12, and 24 for each series
	spi_data <- list()
	spi_data[["scale_3"]] <- SPEI::spi(as.matrix(file_data$data), scale = 3, verbose = FALSE)$fitted
	spi_data[["scale_3"]][is.na(spi_data[["scale_3"]])] = 0
	spi_data[["scale_3"]][is.infinite(spi_data[["scale_3"]]) & spi_data[["scale_3"]] > 0] = max(spi_data[["scale_3"]][!is.infinite(spi_data[["scale_3"]])], na.rm = TRUE)
	spi_data[["scale_3"]][is.infinite(spi_data[["scale_3"]]) & spi_data[["scale_3"]] < 0] = min(spi_data[["scale_3"]][!is.infinite(spi_data[["scale_3"]])], na.rm = TRUE)
	spi_data[["scale_3"]][c(1:2), ] = NA
	spi_data[["scale_12"]] <- SPEI::spi(as.matrix(file_data$data), scale = 12, verbose = FALSE)$fitted
	spi_data[["scale_12"]][is.infinite(spi_data[["scale_12"]]) & spi_data[["scale_12"]] > 0] = max(spi_data[["scale_12"]][!is.infinite(spi_data[["scale_12"]])], na.rm = TRUE)
	spi_data[["scale_12"]][is.infinite(spi_data[["scale_12"]]) & spi_data[["scale_12"]] < 0] = min(spi_data[["scale_12"]][!is.infinite(spi_data[["scale_12"]])], na.rm = TRUE)
	rownames(spi_data[["scale_3"]]) <- rownames(file_data$data)
	rownames(spi_data[["scale_12"]]) <- rownames(file_data$data)

	################ Generate arrays and create trend figures
	mkTrend_pval <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_pval)

	# Corrected version of pval
	mkTrend_pval_grid_points <- list()
	for(name in names(mkTrend_pval)){
  	mkTrend_pval_grid_points[[name]] <- stats::p.adjust(mkTrend_pval[[name]], "BH") # Correction explained in the paper "The Stippling Shows Statistically Significant Grid Points"
	}

	mkTrend_slp <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_slp)

	percentage <- calc_data_year_month_station(data = file_data$data, calc_function = calc_percentage)

	# Averages
	data_mean <- calc_data_year_month_station(data = file_data$data, calc_function = base::mean)

	# Coefficients of variation, annual and seasonal standard deviation
	cv <- calc_data_year_month_station(data = file_data$data, calc_function = coef_var)

	dry_spell_trend_data <- list()
	dry_spell_trend_data[["scale_3"]] <- t(apply(spi_data[["scale_3"]], c(2), dry_spell_trend, threshold = 0))
	colnames(dry_spell_trend_data[["scale_3"]]) <- c("p.dur", "p.mag", "trend.dur", "trend.mag")
	dry_spell_trend_data[["scale_12"]] <- t(apply(spi_data[["scale_12"]][c(13:dim(spi_data[["scale_12"]])[1]), , drop = FALSE], c(2), dry_spell_trend, threshold = 0))
	colnames(dry_spell_trend_data[["scale_12"]]) <- c("p.dur", "p.mag", "trend.dur", "trend.mag")

	if(calc_mobile_trends_data){ # It only calculates it for some dates
		mobile_trends_data <- calc_data_year_month_station(data = file_data$data, calc_function = mobile_trends)
	}else{
		mobile_trends_data <- NA
	}

	# Indicates with a 1 that the original data was not NA and with a 0 that it was NA
  data_change <- data_ori[rownames(file_data$data), colnames(file_data$data), drop = FALSE]
  data_change[!is.na(data_change)] <- 1
  data_change[is.na(data_change)] <- 0

  data_change <- as.data.frame(data_change)
  data_change["year"] = as.numeric(substr(rownames(data_change), 8, 11))
  data_change["month"] = as.numeric(months(chron::chron(rownames(data_change), format=c(dates = "d/m/yy", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s"))))
  data_change = data_change[c("year", "month", colnames(data_change)[1:(dim(data_change)[2]-2)])]

	return(list(coor = file_data$coor, data_change = data_change, regional_series = average, pvalue_trend = mkTrend_pval, pvalue_corrected_trend = mkTrend_pval_grid_points, magnitude_change_mm = mkTrend_slp, magnitude_change_percentage = percentage, average_precipitation_stations = data_mean, cv_precipitation_stations = cv, dry_spell_trend_data = dry_spell_trend_data, mobile_trends_data = mobile_trends_data, change_values_homog = file_data$change_values, break_points_homog = file_data$break_points, pvalue_trend_homog = file_data$mkTrend_pval_pre_homogenize, magnitude_change_mm_homog = file_data$mkTrend_slp_homogenize, magnitude_change_percentage_homog = file_data$percentage_homogenize))
}

#' Final output with all statistics, regional average series, trends, SPI...
#'
#' @param data data and coordinates 
#' @param data_ori original data
#'
#' @return data and coordinates
#' @export
#'
calculate_statistics <- function(data, data_ori){
	suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = length(data)))

	snowfall::sfExport(list = c("read_years", "calc_data_year_month_station", "calc_mkTrend_pval", "calc_data_year", "mkTrend", "calc_mkTrend_slp", "calc_percentage", "coef_var", "dry_spell_trend", "mobile_trends"), local = TRUE, namespace = "mediterraneancalculations")
	# snowfall::sfLibrary("mediterraneancalculations", character.only = TRUE) # Quitar

	# calculate_statistics_data(file_data = data[[1]], data_ori = data_ori)
	data_save <- snowfall::sfLapply(data, calculate_statistics_data, data_ori)
	suppressMessages(snowfall::sfStop())
	return(data_save)
}

#' Below 28 degrees north, remove stations
#'
#' @param data data and coordinates
#'
#' @return data and coordinates
#' @export
#'
delete_zones <- function(data){
	min_north <- 28 # latitude

	i_ini <- names(data)[length(data)]
	for(i_ini in names(data)){
		data[[i_ini]]$coor <- data[[i_ini]]$coor[data[[i_ini]]$coor[, "lat"] >= min_north, , drop = FALSE]
		data[[i_ini]]$data <- data[[i_ini]]$data[, rownames(data[[i_ini]]$coor), drop = FALSE]
	}
	return(data)
}

#' Reads precipitation files, calculates statistics, and saves the results
#' The input files are 2 CSVs: one with coordinates in degrees (stations in rows and lat/lon in columns) and the other with monthly data (dates in rows and stations in columns)
#'
#' @param file_data path to the data file
#' @param file_coor path to the coordinates file
#'
#' @return None
#' @export
#'
main_mediterranean_calculations <- function(file_data, file_coor){
	pb <- utils::txtProgressBar(min = 0,  # Minimum value of the progress bar
                     max = 100,  # Maximum value of the progress bar
                     style = 3,  # Style of the bar (other styles: style = 1 and style = 2)
                     width = 50, # Width of the bar. Default: getOption("width")
                     char = "=")

	folder <- alexanderson_folder

	# Read the data
	read_all_data <- read_data(file_data, file_coor)

	utils::setTxtProgressBar(pb, 5)

	data_statistics <- main_mediterranean_calculations_(read_all_data, folder)

	utils::setTxtProgressBar(pb, 95)

	saveRDS(data_statistics, file = "mediterranean_calculations.rds")	
	# data_save <- readRDS("mediterranean_calculations.rds")

	utils::setTxtProgressBar(pb, 100)

	# save_csvs(i_ini = "prueba", folder_name = "prueba", data_save = data_first_fill$start_1981$data, coor_save = data_first_fill$start_1981$coor)
	# save.image("main_mediterranean_calculations.RData")
	# load("main_mediterranean_calculations.RData")
	return(data_statistics)
}

#' Calculates statistics for a country
#'
#' @param read_all_data input data
#' @param folder directory where files are saved
#' @param pb progress bar
#'
#' @return None
#' @export
#'
main_mediterranean_calculations_ <- function(read_all_data, folder, pb = NULL){
	max_dist_eval <- 100
	max_dist_fill <- 200
	max_dist_second_eval <- 200

	# A) Quality control
	control_data <- mediterranean_calculations(data = read_all_data, max_dist_eval = max_dist_eval)
	statistics_data_removed <- save_delete_data(ori_data = read_all_data$data, process_data = control_data$data, folder = folder)

	if(!is.null(pb)) { utils::setTxtProgressBar(pb, 15) }

	# B) Reconstruction - First fill and second fill, only between the data already selected in each dataset
	fill_data <- fill_series(control_data = control_data, min_correlation = min_correlation, max_dist = max_dist_fill)
	data_first_fill <- save_data(data_ori = read_all_data$data_ori, control_data = fill_data)
	data_second_fill <- second_data_fill(data = data_first_fill, max_dist_eval = max_dist_second_eval)

	# Remove stations below 28 degrees north
	data_second_fill_zone <- delete_zones(data = data_second_fill)

	if(!is.null(pb)) { utils::setTxtProgressBar(pb, 35) }

	# C) Homogenization - Alexanderson test
	data_homogenize <- alexanderson_homogenize(data = data_second_fill_zone, folder = folder)

	if(!is.null(pb)) { utils::setTxtProgressBar(pb, 55) }

	# D) Analysis - Final output: regional average series, trends, and SPI
	data_statistics <- calculate_statistics(data = data_homogenize, data_ori = read_all_data$data_ori)

	if(!is.null(pb)) { utils::setTxtProgressBar(pb, 93) }

	data_statistics$statistics_data_removed <- statistics_data_removed
	data_statistics$availability_data_original <- list()
	data_statistics$availability_data_original$data <- read_all_data$data
	data_statistics$availability_data_original$coor <- read_all_data$coor
	data_statistics$availability_data_original$data[!is.na(data_statistics$availability_data_original$data)] <- 1
	data_statistics$availability_data_original$data[is.na(data_statistics$availability_data_original$data)] <- 0
	rows_availability_data_original <- chron::chron(rownames(data_statistics$availability_data_original$data), format=c(dates = "d/m/yy", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s")) <= chron::chron(paste0("31/12/", i_end), format=c(dates = "d/m/y", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s"))
	data_statistics$availability_data_original$data <- data_statistics$availability_data_original$data[rows_availability_data_original, ]

	return(data_statistics)
}
