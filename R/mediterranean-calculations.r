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


#' Lee los ficheros de precipitación, hace el relleno y guarda los ficheros de salida
#' Los ficheros de entrada son 2 CSVs uno de coordenadas en grados (filas las estaciones y columnas lat y lon y otro de datos mensuales con fechas en filas y las estaciones en las columnas)
#' Control de calidad: Se estaciones con menos de 20 años de datos y usando las 10 más correlacionadas a menos de 200 km, se desechan los datos con un percentil de diferencia de más de 0.6.
#' Rellenado mensual de las series: Usamos estaciones a menos de 200km con correlación por encima de 0.7. Excepto para junio, julio y agosto, meses para los que rellenamos con la más cercana.
#' Homogeneiza, test de Alexandersson (snht)
#' Guarda la salida en 5 fichero con los datos
#'  5 ficheros que indican si cada datos es original o rellenado 
#'  y 5 ficheros de coordenadas para las estaciones de cada fichero de datos, que son:
#' - 1870 a 2020 con más de 80 años originales
#' - 1900 a 2020 con más de 80 años originales
#' - 1930 a 2020 con más de 60 años originales
#' - 1950 a 2020 con más de 40 años originales
#' - 1990 a 2020 con más de 30 años originales
#'
#' @param file_data ruta del fichero de datos
#' @param file_coor ruta del fichero de coordenadas
#' @param max_dist máxima distancia entre 2 estaciones para usar una para evaluar o completar la otra
#'
#' @return data y coor con los datos que no pasan el control eliminados
#' @export
#'
mediterranean_calculations <- function(file_data = "data.csv", file_coor = "coor.csv", max_dist){

	read_all_data <- read_data(file_data, file_coor)

	####################

	# Control de calidad, se eliminan series completas y datos sueltos inválidos
	control_data <- quality_control(data = read_all_data$data, coor = read_all_data$coor, max_dist = max_dist)
	data_select = read_all_data$data[rownames(control_data$data), colnames(control_data$data)]

	#################### save.image("quality_control.RData")

	# Rellenado de series
	fill_data <- fill_series(control_data, min_correlation, max_dist)

	####################

	# Guardar los datos
	data <- save_data(data_ori = read_all_data$data_ori, control_data = fill_data)
	# save_parameters(folder_name, max_dist)
	# save.image("calculos_mediterraneo.RData") #fergus: para pruebas
	return(data)
}

#' Realiza un segundo relleno
#' Para cada base de datos (1870, 1900...)
#' Estaciones con más del 90 o 95% de datos rellenos
#' Ordenamos las estaciones por correlación (mínima 0.5) y rellenamos usando los 10 métodos...
#' Las estaciones sin relleno total las tiramos
#'
#' @param file_data ruta del fichero de datos
#' @param fillable_years años rellenables con la media mensual de la propia estación
#'
#' @return None
#' @export
#'
second_data_fill_data <- function(file_data, fillable_years = 36){

	# Retiramos estaciones sin 95% de datos
	no_nas <- dim(file_data$data)[1] - apply(file_data$data, c(2), sum_no_nas) < dim(file_data$data)[1] * percentage_filled_data / 100
	control_data <- list(data = file_data$data[, no_nas], coor = file_data$coor[no_nas, ])

	if(sum(no_nas) > 0){
		# segundo relleno
		fill_data <- fill_series(control_data = control_data, min_correlation = min_second_correlation, max_dist = NA)

		# Retiramo estaciones no completas
		no_nas <- dim(fill_data$data)[1] - apply(fill_data$data, c(2), sum_no_nas) == 0
		all_data <- list(data = fill_data$data[, no_nas], coor = fill_data$coor[no_nas, ])
		if(length(all_data$data) > 0 & is.null(dim(all_data$data))){
			all_data$data = array(all_data$data, dim = c(length(all_data$data), 1), dimnames = list(rownames(fill_data$data), colnames(fill_data$data)[no_nas]))
		}

		# Si no tenemos datos, intentamos rellenar alguna serie consigo misma
		if(dim(all_data$coor)[1] <= 6){
			# save.image("second_data_fill_data_fill_series.RData") # fergus: pruebas
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
#'
#' @return data y coor
#' @export
#'
second_data_fill <- function(data){
  all_data <- list()
  i_ini <- i_inis[length(i_inis)]
  for(i_ini in i_inis){
  	file_data <- data[[paste0("start_", i_ini)]]

  	# Años rellenables con la media mensual de la propia estación
  	if(i_ini == 1980){
  		fillable_years <- 24 # 2 años
  	}else{
  		fillable_years <- 36 # 3 años
  	}

  	all_data[[paste0("start_", i_ini)]] <- second_data_fill_data(file_data = file_data, fillable_years = fillable_years)  
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
#'
#' @param file_data ruta del fichero de datos
#'
#' @return None
#' @export
#'
alexanderson_homogenize_data <- function(file_data){
	data_save <- file_data$data
	all_series <- colnames(file_data$data)

	# Diferencias de un día respecto al anterior
	differences_0 <- file_data$data
	differences_1 <- rbind(rep(0, dim(differences_0)[2]), differences_0[c(1:(dim(differences_0)[1]-1)), ])
	differences <- differences_0 - differences_1

 	data_cor_list <- near_correlations(data = differences, coor = file_data$coor, max_dist = NA)
	months <- names(data_cor_list)
	i_month <- months[2]
	for(i_month in months){
	    data_cor <- data_cor_list[[i_month]]
    	date_month <- grepl(i_month, rownames(file_data$data))
    	data_save_month <- file_data$data[date_month, ]
    	data_month <- NA
    	i_snht <- 1
		while(i_snht < 20 & (length(data_month) == 1 | sum(data_save_month != data_month, na.rm = TRUE) > 0)){
			# print(paste(i_month))
			data_month <- data_save_month
			i_series <- all_series[1]
			for(i_series in all_series){
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
			    		q_snht <- snht(cbind(1:length(Q), Q), 95) # https://sites.google.com/site/phaenggi23/rscripts
			    		break_serie <- q_snht$T0x
	      				pre_post_mean <- round(mean(data_month[c((break_serie + 1):(dim(data_month)[1])), i_series]) / mean(data_month[c(1:break_serie), i_series]), digits = 15)
	      				pre_post_mean[pre_post_mean == Inf] = NA
			    		pre_post_mean[pre_post_mean == -Inf] = NA
			    		pre_post_mean[is.nan(pre_post_mean)] = NA
			    		if (q_snht$T0 > q_snht$Tc && !is.na(pre_post_mean) && pre_post_mean != 1 && pre_post_mean != 0) {       			
		        			data_save_month[c(1:break_serie), i_series] <- pre_post_mean * data_month[c(1:break_serie), i_series]
	    				} #if
					} #if
				} #if
			} #for
				i_snht <- i_snht + 1
		} #while
		data_save[date_month, ] <- data_save_month
	}#for
	return(data_save)
}

#' Test de Alexanderson para todos los ficheros disponibles (que han pasado el segundo relleno con éxito)
#'
#' @param data data y coor 
#'
#' @return data y coor
#' @export
#'
alexanderson_homogenize <- function(data){
	dir.create(alexanderson_folder, recursive = TRUE, showWarnings = FALSE)
  
  data_save <- list()
	i_ini <- i_inis[length(i_inis)]
	for(i_ini in i_inis){
    	file_data <- data[[paste0("start_", i_ini)]]
    	if(length(file_data) > 1 && dim(file_data$coor)[1] > 6){
	    	data_save[[paste0("start_", i_ini)]] <- list(data = alexanderson_homogenize_data(file_data), coor = file_data$coor)
				save_csvs(i_ini, folder_name = alexanderson_folder, data_save = data_save[[paste0("start_", i_ini)]]$data, coor_save = data_save[[paste0("start_", i_ini)]]$coor)
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
#' @param file_data ruta del fichero de datos
#' @param calc_mobile_trends_data Calculamos mobile_trends_data o no
#'
#' @return None
#' @export
#'
calculate_statistics_data <- function(file_data, calc_mobile_trends_data = TRUE){
	# Serie promedio de todo el país
	average <- apply(file_data$data, c(1), mean)

	# SPI a escalas 3, 9, 12, y 24 de cada serie
	spi_data <- list()
	spi_data[["scale_3"]] <- SPEI::spi(as.matrix(file_data$data), scale = 3, verbose = FALSE)$fitted
	# spi_data[["scale_3"]][is.infinite(spi_data[["scale_3"]]) | is.na(spi_data[["scale_3"]])] = 0
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

	mkTrend_slp <- calc_data_year_month_station(data = file_data$data, calc_function = calc_mkTrend_slp)

	percentage <- calc_data_year_month_station(data = file_data$data, calc_function = calc_percentage)

	dry_spell_trend_data <- list() # fergus:scale 3 tiene NAs porque SPI 3 tiene NAs????
	dry_spell_trend_data[["scale_3"]] <- t(apply(spi_data[["scale_3"]], c(2), dry_spell_trend, threshold = 0))
	colnames(dry_spell_trend_data[["scale_3"]]) <- c("p.dur", "p.mag", "trend.dur", "trend.mag")
	dry_spell_trend_data[["scale_12"]] <- t(apply(spi_data[["scale_12"]][c(13:dim(spi_data[["scale_12"]])[1]), ], c(2), dry_spell_trend, threshold = 0))
	colnames(dry_spell_trend_data[["scale_12"]]) <- c("p.dur", "p.mag", "trend.dur", "trend.mag")

	if(calc_mobile_trends_data){ #Solo lo calculamos para algunas fechas
		mobile_trends_data <- calc_data_year_month_station(data = file_data$data, calc_function = mobile_trends)
	}else{
		mobile_trends_data <- NA
	}

	return(list(coor = file_data$coor, average = average, pval = mkTrend_pval, trend = mkTrend_slp, percentage = percentage, dry_spell_trend_data = dry_spell_trend_data, mobile_trends_data = mobile_trends_data))
}

#' Salida final: serie regional promedio, tendencias, SPI...
#' Un RData o RDS y ficheros CSV sueltos por base da datos
#' Identificador, X, Y, series...
#'
#' @param data data y coor 
#'
#' @return data and coor
#' @export
#'
calculate_statistics <- function(data){
	data_save <- list()
	i_ini <- i_inis[length(i_inis)]
	for(i_ini in i_inis){
    	file_data <- data[[paste0("start_", i_ini)]]
  		if(i_ini <= 1900){ # Calcular mobile_trends_data es muy lento, lo haremos solo en algunos casos
    		calc_mobile_trends_data = TRUE
			}else{
  			calc_mobile_trends_data = FALSE
  		}
    	if(sum(!is.na(file_data)) > 0){
	    	data_save[[paste0("start_", i_ini)]] <- calculate_statistics_data(file_data, calc_mobile_trends_data)
    	}
  	}
	saveRDS(data_save, file = "mediterranean_calculations.rds")	
	# data_save <- readRDS("mediterranean_calculations.rds")
	return(data_save)
}

#' Lee los ficheros de precipitación, calcula estadísitocs y guarda los resultados
#' Los ficheros de entrada son 2 CSVs uno de coordenadas en grados (filas las estaciones y columnas lat y lon y otro de datos mensuales con fechas en filas y las estaciones en las columnas)
#'
#' @param file_data ruta del fichero de datos
#' @param file_coor ruta del fichero de coordenadas
#' @param max_dist máxima distancia entre 2 estaciones para usar una para evaluar o completar la otra
#'
#' @return None
#' @export
#'
main_mediterranean_calculations <- function(file_data, file_coor, max_dist){
	# Control de calidad y primer relleno
	data_first_fill <- mediterranean_calculations(file_data, file_coor, max_dist)

	# Relleno solo entre los datos ya seleccionados en cada base de datos
	data_second_fill <- second_data_fill(data = data_first_fill)

	# Homogeneizar - test de Alexanderson
	data_homogenize <- alexanderson_homogenize(data = data_second_fill)

	# Salida final: serie regional promedio, tendencias y SPI
	data_statistics <- calculate_statistics(data = data_homogenize)

	# save_csvs(i_ini = "prueba", folder_name = "prueba", data_save = data_first_fill$start_1980$data, coor_save = data_first_fill$start_1980$coor)
	# save.image("main_mediterranean_calculations.RData")
	# load("main_mediterranean_calculations.RData")
	return(data_statistics)
}
