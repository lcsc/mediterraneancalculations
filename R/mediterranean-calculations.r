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


#' Hace el control de calidad
#' Control de calidad: Se estaciones con menos de 20 años de datos y usando las 10 más correlacionadas a menos de 200 km, se desechan los datos con un percentil de diferencia de más de 0.6.
#'
#' @param data ruta del fichero de datos
#' @param max_dist_eval máxima distancia entre 2 estaciones para usar una para evaluar una con la otra
#'
#' @return data y coor con los datos que no pasan el control eliminados
#' @export
#'
mediterranean_calculations <- function(data, max_dist_eval){

	# Eliminamos datos si tenemos 5 meses o más seguidos de 0s, si uno de los meses implicados tiene menos del 70% de ceros
	delete_zero_data <- delete_zero(data = data$data)

	# Control de calidad, se eliminan series completas y datos sueltos inválidos
	control_data <- quality_control(data = delete_zero_data, coor = data$coor, max_dist = max_dist_eval, max_diff_anomaly = max_diff_anomaly, max_diff_anomaly_0 = max_diff_anomaly_0)

	return(control_data)
}

#' Realiza un segundo relleno
#' Para cada base de datos (1870, 1900...)
#' Estaciones con más del 90 o 95% de datos rellenos
#' Ordenamos las estaciones por correlación (mínima 0.5) y rellenamos usando los 10 métodos...
#' Las estaciones sin relleno total las tiramos
#'
#' @param file_data ruta del fichero de datos
#' @param fillable_years años rellenables con la media mensual de la propia estación
#' @param max_dist máxima distancia permitida para el relleno
#'
#' @return None
#' @export
#'
second_data_fill_data <- function(file_data, fillable_years = 36, max_dist = NA){

	# Retiramos estaciones sin 95% de datos
	no_nas <- (dim(file_data$data)[1] - apply(file_data$data, c(2), sum_no_nas)) < (dim(file_data$data)[1] * (100 - percentage_filled_data) / 100)
	control_data <- list(data = file_data$data[, no_nas], coor = file_data$coor[no_nas, ])

	if(sum(no_nas) > 0){
		# segundo relleno
		fill_data <- fill_series(control_data = control_data, min_correlation = min_second_correlation, max_dist = max_dist)

		# Retiramo estaciones no completas
		no_nas <- dim(fill_data$data)[1] - apply(fill_data$data, c(2), sum_no_nas) == 0
		all_data <- list(data = fill_data$data[, no_nas], coor = fill_data$coor[no_nas, ])
		if(length(all_data$data) > 0 & is.null(dim(all_data$data))){
			all_data$data = array(all_data$data, dim = c(length(all_data$data), 1), dimnames = list(rownames(fill_data$data), colnames(fill_data$data)[no_nas]))
		}

		# Si no tenemos datos, intentamos rellenar alguna serie consigo misma
		if(dim(all_data$coor)[1] <= 6){
			# Rellenar 2 o 3 años de datos con la media mensual de la propia estación
			fill_data <- fill_unfillable_station(data = fill_data, fillable_years = fillable_years)

			# Retiramo estaciones no completas
			no_nas <- dim(fill_data$data)[1] - apply(fill_data$data, c(2), sum_no_nas) == 0
			all_data <- list(data = fill_data$data[, no_nas], coor = fill_data$coor[no_nas, ])
			if(length(all_data$data) > 0 & is.null(dim(all_data$data))){
				all_data$data = array(all_data$data, dim = c(length(all_data$data), 1), dimnames = list(rownames(fill_data$data), colnames(fill_data$data)[no_nas]))
			}
		}

	}else{
		all_data <- NA
	}
	return(all_data)
}

#' Realiza un segundo relleno
#' Para cada base de datos (1870, 1900...)
#' Estaciones con más del 90 o 95% de datos rellenos
#' Ordenamos las estaciones por correlación (mínima 0.5) y rellenamos usando los 10 métodos...
#' Las estaciones sin relleno total las tiramos
#'
#' @param data data y coor 
#' @param max_dist_eval maxima distancia para el relleno
#'
#' @return data y coor
#' @export
#'
second_data_fill <- function(data, max_dist_eval = NA){
  all_data <- list()
  i_ini <- names(data)[length(data)]
  for(i_ini in names(data)){
  	file_data <- data[[i_ini]]

  	# Años rellenables con la media mensual de la propia estación
  	if(i_ini == paste0("start_", i_inis[length(i_inis)])){
  		fillable_years <- 24 # 2 años
  	}else{
  		fillable_years <- 36 # 3 años
  	}

  	data_fill <- second_data_fill_data(file_data = file_data, fillable_years = fillable_years, max_dist = max_dist_eval)
  	if(sum(!is.na(data_fill)) > 0){
  		all_data[[i_ini]] <- data_fill
  	}
	}
	return(all_data)
}

#' Homogeneizar - test de Alexanderson
#' Lo usamos en code_web_maps/snht_functions.R
#' Existe también librería snht de R
#' -- Serie de referencia, compara y corrige
#' Para cada base de datos (1870, 1900...)
#' Elegimos las 5 series más correlacionadas usandos las serie de diferencias
#' Con las 5 hacemos una media ponderada, (correlación * dato1 + ...) / sum(correlaciones) y será la serie de referencia
#' Alexanderson nos dará un punto de ruptura y un valor ratio por el que multiplicar la parte antigua... iterar mientras de puntos de ruptura
#' Guardar estadísticos de inhomogeneidades. Básicamente número de datos cambiados en cada series y momento de la inhomogeneidad. 
#' CSV con número de datos cambiados y CSV con punto de inhomogeneidad - todo x 12 meses
#'
#' @param file_data ruta del fichero de datos
#' @param no_use_series series que no se homogeneizarán
#'
#' @return None
#' @export
#'
alexanderson_homogenize_data <- function(file_data, no_use_series = c()){

	max_i_snht <- 20
	significance_level <- 100 * (1 - 0.01)
	data_save <- file_data$data

	# Buscamos inhomogenidades solo en las series que no estaban ya en los periodos previos
	all_series <- colnames(file_data$data)
	use_series <- all_series[!all_series %in% no_use_series]

	# Diferencias de un día respecto al anterior
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

		    		# Media ponderada
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

#' Test de Alexanderson para todos los ficheros disponibles (que han pasado el segundo relleno con éxito)
#'
#' @param data data y coor
#' @param folder directorio para guardar los datos de salida
#'
#' @return data y coor
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
    	if(length(file_data) > 1 && dim(file_data$coor)[1] > 6){
    		# Las estaciones que ya estaban en el grup precio han de ser exáctamente las mismas
    		if(!is.null(file_data_pre)){
    			delete_stations <- pre_stations[!pre_stations %in% colnames(file_data$data)]
    			if(length(delete_stations) > 0){
    				file_data$data <- cbind(file_data$data, file_data_pre$data[rownames(file_data$data), delete_stations])
    				file_data$coor <- rbind(file_data$coor, file_data_pre$coor[delete_stations, ])
    			}
    			file_data$data[, pre_stations] <- file_data_pre$data[rownames(file_data$data), pre_stations]
    		}
    		homogenize_data <- alexanderson_homogenize_data(file_data = file_data, no_use_series = pre_stations)
	    	
    		# Generamos estadísticos pre homogenización
    		mkTrend_pval <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_pval)

  			mkTrend_slp <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_slp)

				percentage <- calc_data_year_month_station(data = file_data$data, calc_function = calc_percentage)

	    	data_save[[i_ini]] <- list(data = homogenize_data$data, coor = file_data$coor, change_values = homogenize_data$change_values, break_points = homogenize_data$break_points, mkTrend_pval_pre_homogenize = mkTrend_pval, mkTrend_slp_homogenize = mkTrend_slp, percentage_homogenize = percentage)

	    	pre_stations <- rownames(data_save[[i_ini]]$coor)
	    	file_data_pre <- data_save[[i_ini]]
				save_csvs(i_ini, folder_name = folder, data_save = data_save[[i_ini]]$data, coor_save = data_save[[i_ini]]$coor)
  		}#if
	}#for
	return(data_save)
}

#' Calcula estadísticos de los datos
#' Tendencia mensual, estacional y anual, paquete Trend, función sens.slope
#' Sumar 1 a todo para evitar 0s
#' Significación, paquete modifiedmk, funsión bbsmk
#' Serie promedio de todo el país
#' SPI a escalas 3, 12, y 24 de cada serie, importante que sea imposible invertir las operaciones
#' Código Sergio para generar arrays y hacer figuras de tendencia
#'
#' @param file_data datos y coordenadas
#' @param data_ori datos originales
#'
#' @return None
#' @export
#'
calculate_statistics_data <- function(file_data, data_ori){

	mobile_trends_calc <- c(1871, 1901, 1931)
	i_ini <- read_years(rownames(file_data$data))[1]
	if(i_ini %in% mobile_trends_calc){ # Calcular mobile_trends_data es muy lento, lo haremos solo en algunos casos
		calc_mobile_trends_data = TRUE
	}else{
		calc_mobile_trends_data = FALSE
	}

	# Serie promedio de todo el país
	average <- apply(file_data$data, c(1), mean)

	# SPI a escalas 3, 9, 12, y 24 de cada serie
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

	################ Código Sergio para generar arrays y hacer figuras de tendencia
	mkTrend_pval <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_pval)

	# Versión corregida de pval
	mkTrend_pval_grid_points <- list()
	for(name in names(mkTrend_pval)){
  	mkTrend_pval_grid_points[[name]] <- stats::p.adjust(mkTrend_pval[[name]], "BH") # Corección explicada en el paper "The Stippling Shows Statistically Significant Grid Points"
	}

	mkTrend_slp <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_slp)

	percentage <- calc_data_year_month_station(data = file_data$data, calc_function = calc_percentage)

	# Promedios
	data_mean <- calc_data_year_month_station(data = file_data$data, calc_function = base::mean)

	# Coeficientes de variación, desviación estándar anual y estacional
	cv <- calc_data_year_month_station(data = file_data$data, calc_function = coef_var)

	dry_spell_trend_data <- list()
	dry_spell_trend_data[["scale_3"]] <- t(apply(spi_data[["scale_3"]], c(2), dry_spell_trend, threshold = 0))
	colnames(dry_spell_trend_data[["scale_3"]]) <- c("p.dur", "p.mag", "trend.dur", "trend.mag")
	dry_spell_trend_data[["scale_12"]] <- t(apply(spi_data[["scale_12"]][c(13:dim(spi_data[["scale_12"]])[1]), ], c(2), dry_spell_trend, threshold = 0))
	colnames(dry_spell_trend_data[["scale_12"]]) <- c("p.dur", "p.mag", "trend.dur", "trend.mag")

	if(calc_mobile_trends_data){ #Solo lo calculamos para algunas fechas
		mobile_trends_data <- calc_data_year_month_station(data = file_data$data, calc_function = mobile_trends)
	}else{
		mobile_trends_data <- NA
	}

	# Indicamos con un 1 que el dato original era un NA y con un 0 que no era un NA
  data_change <- data_ori[rownames(file_data$data), colnames(file_data$data)]
  data_change[!is.na(data_change)] <- 1
  data_change[is.na(data_change)] <- 0

  data_change <- as.data.frame(data_change)
  data_change["year"] = as.numeric(substr(rownames(data_change), 8, 11))
  data_change["month"] = as.numeric(months(chron::chron(rownames(data_change), format=c(dates = "d/m/yy", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s"))))
  data_change = data_change[c("year", "month", colnames(data_change)[1:(dim(data_change)[2]-2)])]

	return(list(coor = file_data$coor, data_change = data_change, regional_series = average, pvalue_trend = mkTrend_pval, pvalue_corrected_trend = mkTrend_pval_grid_points, magnitude_change_mm = mkTrend_slp, magnitude_change_percentage = percentage, average_precipitation_stations = data_mean, cv_precipitation_stations = cv, dry_spell_trend_data = dry_spell_trend_data, mobile_trends_data = mobile_trends_data, change_values_homog = file_data$change_values, break_points_homog = file_data$break_points, pvalue_trend_homog = file_data$mkTrend_pval_pre_homogenize, magnitude_change_mm_homog = file_data$mkTrend_slp_homogenize, magnitude_change_percentage_homog = file_data$percentage_homogenize))
}

#' Salida final con todos los estadíticos, serie regional promedio, tendencias, SPI...
#'
#' @param data data y coor 
#' @param data_ori data original
#'
#' @return data and coor
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

#' Por debajo de 28 grados norte, eliminar estaciones
#'
#' @param data data y coor 
#'
#' @return data and coor
#' @export
#'
delete_zones <- function(data){
	min_north <- 28 # latitud

	i_ini <- names(data)[length(data)]
	for(i_ini in names(data)){
		data[[i_ini]]$coor <- data[[i_ini]]$coor[data[[i_ini]]$coor[, "lat"] >= min_north, ]
		data[[i_ini]]$data <- data[[i_ini]]$data[, rownames(data[[i_ini]]$coor)]
	}
	return(data)
}

#' Lee los ficheros de precipitación, calcula estadísticos y guarda los resultados
#' Los ficheros de entrada son 2 CSVs uno de coordenadas en grados (filas las estaciones y columnas lat y lon y otro de datos mensuales con fechas en filas y las estaciones en las columnas)
#'
#' @param file_data ruta del fichero de datos
#' @param file_coor ruta del fichero de coordenadas
#'
#' @return None
#' @export
#'
main_mediterranean_calculations <- function(file_data, file_coor){
	max_dist_eval <- 100
	max_dist_fill <- 200
	max_dist_second_eval <- 200

	# Leemos los datos
	read_all_data <- read_data(file_data, file_coor)

	# A) Control de calidad
	control_data <- mediterranean_calculations(data = read_all_data, max_dist_eval = max_dist_eval)
	statistics_data_removed <- save_delete_data(ori_data = read_all_data$data, process_data = control_data$data, folder = alexanderson_folder)

	# B) Reconstrucción - Primer relleno y segundo relleno solo entre los datos ya seleccionados en cada base de datos
	fill_data <- fill_series(control_data = control_data, min_correlation = min_correlation, max_dist = max_dist_fill)
	data_first_fill <- save_data(data_ori = read_all_data$data_ori, control_data = fill_data)
	data_second_fill <- second_data_fill(data = data_first_fill, max_dist_eval = max_dist_second_eval)

	# Calcular estadísticos de la reconstrucción
	# reconstruction_statistics <- calculate_reconstruction_statistics(sim = data_second_fill, obs = read_all_data$data_ori)

	# Eliminar series por encima debajo de 28 grados Norte
	data_zone_second_fill <- delete_zones(data = data_second_fill)

	# C) Homogeneización - test de Alexanderson
	data_homogenize <- alexanderson_homogenize(data = data_zone_second_fill, folder = alexanderson_folder)

	# D) Análisis - Salida final: serie regional promedio, tendencias y SPI
	data_statistics <- calculate_statistics(data = data_homogenize, data_ori = read_all_data$data_ori)

	data_statistics$statistics_data_removed <- statistics_data_removed
	data_statistics$availability_data_original <- list()
	data_statistics$availability_data_original$data <- read_all_data$data
	data_statistics$availability_data_original$coor <- read_all_data$coor
	data_statistics$availability_data_original$data[!is.na(data_statistics$availability_data_original$data)] <- 1
	data_statistics$availability_data_original$data[is.na(data_statistics$availability_data_original$data)] <- 0
	rows_availability_data_original <- chron::chron(rownames(data_statistics$availability_data_original$data), format=c(dates = "d/m/yy", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s")) <= chron::chron(paste0("31/12/", i_end), format=c(dates = "d/m/y", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s"))
	data_statistics$availability_data_original$data <- data_statistics$availability_data_original$data[rows_availability_data_original, ]

	saveRDS(data_statistics, file = "mediterranean_calculations.rds")	
	# data_save <- readRDS("mediterranean_calculations.rds")

	# save_csvs(i_ini = "prueba", folder_name = "prueba", data_save = data_first_fill$start_1981$data, coor_save = data_first_fill$start_1981$coor)
	# save.image("main_mediterranean_calculations.RData")
	# load("main_mediterranean_calculations.RData")
	return(data_statistics)
}
