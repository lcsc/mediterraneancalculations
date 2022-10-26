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

#' @importFrom stats cor ecdf qnorm shapiro.test ecdf lm coef quantile acf pnorm median na.omit sd
#' @importFrom utils read.table write.table
NULL

#' @importFrom chron chron
#' @importFrom sp CRS coordinates proj4string spTransform
#' @importFrom SpatialTools dist1
#' @importFrom lmom samlmu pelexp cdfexp pelgam cdfgam pelgev cdfgev pelglo cdfglo pelgpa cdfgpa pelgno cdfgno pelln3 cdfln3 pelnor cdfnor pelpe3 cdfpe3 pelwei cdfwei
#' @importFrom SPEI spi
NULL

### Variables globales
crs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #https://spatialreference.org/ref/epsg/wgs-84/ "epsg:4326"
crs84m <- sp::CRS("+proj=geocent +datum=WGS84 +units=m +no_defs") #https://epsg.io/4978 "EPSG:4978"

min_years <- 10 # Años mínimos de datos para que no se elimine una estación
min_overlap <- 10 # Años mínimos de datos de solape entre 2 series para rellenar una con la otra
n_reference_stations <- 10 # Número de datos usados para sacar la serie de referencia en el control de calidad
min_correlation <- 0.7 #Correlación mínima para usar el dato en el relleno
max_diff_anomaly <- 0.6 #Máxima diferencia de anomalías para mantener dato en el control
max_diff_anomaly_0 <- 0.5 #Máxima diferencia de anomalías para mantener dato en el control, si el dato es 0
percentage_filled_data <- 95 #Número de datos que tiene una estación para entrar en el segundo relleno
min_second_correlation <- 0.5 #Correlación mínima entre estaciones para usarse en el segundo relleno

####################

# library(rworldmap, include.only = "getMap")
# world = rworldmap::getMap(resolution = "high")

time_format <- c(dates = "d/m/yy", times = "h:m:s")

pel_functions <- list(exp = lmom::pelexp, gam = lmom::pelgam, gev = lmom::pelgev, glo = lmom::pelglo, gpa = lmom::pelgpa, gno = lmom::pelgno, ln3 = lmom::pelln3, nor = lmom::pelnor, pe3 = lmom::pelpe3, wei = lmom::pelwei)
cdf_functions <- list(exp = lmom::cdfexp, gam = lmom::cdfgam, gev = lmom::cdfgev, glo = lmom::cdfglo, gpa = lmom::cdfgpa, gno = lmom::cdfgno, ln3 = lmom::cdfln3, nor = lmom::cdfnor, pe3 = lmom::cdfpe3, wei = lmom::cdfwei)

alexanderson_folder <- "results"

i_inis <- c(1870, 1880, 1890, 1900, 1930, 1950, 1980)
i_end <- 2020
i_years <- array(c(80, 80, 80, 80, 60, 40, 30), length(i_inis), dimnames = list(i_inis))

####################

#' Ordena los datos y devuelve una lista con el orden
#'
#' @param data datos
#'
#' @return list
#' @export
#'
order_data <- function(data){
  data_order <- apply(data, c(1), order, na.last = NA)
  if(class(data_order)[1] == "matrix"){
    data_order <- as.list(as.data.frame(data_order))
  }
  return(data_order)
}

#' Número de datos distintos de NA
#'
#' @param data datos
#'
#' @return número de datos no NAs
#' @export
#'
sum_no_nas <- function(data){
  return(sum(!is.na(data)))
}

#' Primeros datos válidos (no NAs)
#'
#' @param data datos
#' @param n_reference_stations numero de datos a devolver
#'
#' @return primeros datos distintos de NA
#' @export
#'
select_data <- function(data, n_reference_stations){
  return(data[!is.na(data)][1:n_reference_stations])
}

#' Anomalías de los datos
#'
#' @param data datos mensuales
#'
#' @return anomalías de los datos
#' @export
#'
apply_ecdf_month <- function(data){
  l_mom <- lmom::samlmu(data, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)

  # Seleccionamos la mejor de las 10 funciones
  p.value = array(NA, dim = length(pel_functions), dimnames = list(names(pel_functions)))
  best_funcion = names(pel_functions)[3]
  for(best_funcion in names(pel_functions)){
    shapiro_value <- tryCatch({
      par_exp <- pel_functions[[best_funcion]](l_mom)
      cdf_exp <- cdf_functions[[best_funcion]](data, para = par_exp)
      norm_exp <- stats::qnorm(cdf_exp)
      shapiro_value <- stats::shapiro.test(norm_exp)$p.value
    }, error = function(cond) {
      NA
    })
    p.value[best_funcion] <- shapiro_value
  }
  if(sum(!is.na(p.value)) > 0){
    best_funcion <- names(pel_functions)[which(p.value == max(p.value, na.rm = TRUE))[1]]

    # Utlizamos la mejor de las 10 funciones
    par_exp <- pel_functions[[best_funcion]](l_mom)
    cdf_exp <- cdf_functions[[best_funcion]](data, para = par_exp)
  }else{
    cdf_exp <- rep(NA, length(data))
  }

  return(cdf_exp)
}

#' Anomalías de los datos mensuales
#'
#' @param data datos mensuales
#'
#' @return anomalías de los datos
#' @export
#'
apply_ecdf <- function(data){
  rownames_data <- rownames(data)
  i <- 1
  for(i in c(0:11)){
    data[c(1:length(data))%%12 == i] <- apply_ecdf_month(data[c(1:length(data))%%12 == i])
  }
  rownames(data) <- rownames_data
  return(data)
}

#' Devuelve la correlación entre las estaciones
#' Sin tener en cuenta las que están a más de 200 km (NA en esos casos)
#'
#' @param data datos mensuales
#' @param coor coordenadas de las estaciones que corresponden con data
#' @param max_dist distancia máxima entre las series a utilizar
#'
#' @return correlación entre las estaciones
#' @export
#'
near_correlations <- function(data, coor, max_dist){
  # Pasar coordenadas de grados a metros
  points_data <- as.data.frame(coor)
  sp::coordinates(points_data) <- c('lat', 'lon')
  sp::proj4string(points_data) <- crs84
  # plot(world)
  # plot(points_data, add = TRUE)
  # plot(points_data)
  # plot(world, add = TRUE)
  points_data <- sp::spTransform(points_data, crs84m)

  # Eliminar estaciones a más de 200 km para los cálculos del control de calidad
  dist <- SpatialTools::dist1(sp::coordinates(points_data))
  colnames(dist) <- rownames(dist) <- rownames(coor)
  if(!is.na(max_dist)){
    dist[dist > max_dist * 1000] <- NA
  }
  dist[dist == 0] <- NA

  data_cor_list <- list()
  # Calcular correlaciones (por meses)
  while(length(data_cor_list) <= 11){
    data_month <- data[c(1:dim(data)[1])%%12==length(data_cor_list), ]
    name_month <- substr(rownames(data_month)[1], 4, 6)
    data_cor <- stats::cor(data_month, use = 'pairwise.complete.obs')
    data_cor[is.na(dist)] <- NA
    data_cor_list[[name_month]] <- data_cor
  }
  return(data_cor_list)
}

#' Devuelve las estaciones por orden de cercanía
#' Sin tener en cuenta las que están a más de 200 km (NA en esos casos)
#'
#' @param data datos mensuales
#' @param coor coordenadas de las estaciones que corresponden con data
#' @param max_dist distancia máxima entre las series a utilizar
#'
#' @return correlación entre las estaciones
#' @export
#'
near_estations <- function(data, coor, max_dist){
  # Pasar coordenadas de grados a metros
  points_data <- as.data.frame(coor)
  sp::coordinates(points_data) <- c('lat', 'lon')
  sp::proj4string(points_data) <- crs84
  # plot(world)
  # plot(points_data, add = TRUE)
  # plot(points_data)
  # plot(world, add = TRUE)
  points_data <- sp::spTransform(points_data, crs84m)

  # Eliminar estaciones a más de 200 km para los cálculos del control de calidad
  dist <- SpatialTools::dist1(sp::coordinates(points_data))
  colnames(dist) <- rownames(dist) <- rownames(coor)
  if(!is.na(max_dist)){
    dist[dist > max_dist * 1000] <- NA
  }
  dist[dist == 0] <- NA

  #Ordenar según distancias 
  data_order <- order_data(data = dist)

  return(data_order)
}

#' Control de calidad
#' Estaciones con menos de 20 años de datos retirar
#' Usando las 10 más correlacionadas a menos de 200 km, desechar si promedio de percentil se diferencia en más de 0.6
#'
#' @param data datos  
#' @param coor coordenadas
#' @param max_dist máxima distancia entre 2 estaciones para usar una para evaluar o completar la otra
#'
#' @return data y coor con los datos que no pasan el control eliminados
#' @export
#'
quality_control <- function(data, coor, max_dist){

  # Eliminar estaciones con menos de min_years años de datos
  min_n_data <- min_years * 12
  no_nas <- apply(data, c(2), sum_no_nas) > min_n_data
  data <- data[, no_nas]
  coor <- coor[colnames(data), ]

  # Eliminar estaciones con menos de 1 año de datos para algún mes
  no_nas <- NA
  i <- 1
  for(i in c(0:11)){
    month_no_nas <- apply(data[c(1:dim(data)[1])%%12 == i, ], c(2), sum_no_nas) > 1
    if(length(no_nas) > 1){
      no_nas <- no_nas & month_no_nas
    }else{
      no_nas <- month_no_nas
    }
  }
  data <- data[, no_nas]
  coor <- coor[colnames(data), ]

  all_series <- colnames(data)
  data_correct <- as.matrix(data) 

  #Calcular anomalías de cada dato dentro de cada estación
  data_percent <- apply(data, c(2), apply_ecdf)

  ### Correlación por meses
  data_near <- near_estations(data, coor, max_dist = max_dist)

  months <- unique(substr(rownames(data), 4, 6))
  # Para cada mes
  i_month <- months[1]
  for(i_month in months){
    # Seleccionamos los datos de mes
    data_percent_month <- data_percent[grepl(i_month, rownames(data_percent)), ]

    #Ordenar según distancias
    data_order <- data_near
    if(length(data_order) > 0){
      # 10 series más cercanas
      i_serie <- all_series[1]
      for(i_serie in all_series){
        i_order <- data_order[[i_serie]] #Orden de cercanía con las otras series
        data_ori <- data[rownames(data_percent_month), i_serie]
        data_10 <- data_percent_month[, i_order] # Datos en orden de las otras series
        data_refe_10 <- t(apply(data_10, c(1), select_data, n_reference_stations = n_reference_stations))
        data_refe <- apply(data_refe_10, c(1), mean, na.rm = TRUE)
        data_correct[rownames(data_percent_month)[!is.na(data_ori) & data_ori!=0 & !is.na(data_percent_month[, i_serie]) & !is.na(data_refe) & abs(data_percent_month[, i_serie] - data_refe) > max_diff_anomaly], i_serie] <- NA
        data_correct[rownames(data_percent_month)[!is.na(data_ori) & data_ori==0 & !is.na(data_percent_month[, i_serie]) & !is.na(data_refe) & abs(data_percent_month[, i_serie] - data_refe) > max_diff_anomaly_0], i_serie] <- NA
      }
    }
  }
  return(list(data = data_correct, coor = coor))
}

#' Calcula el tiempo de solape existente entre cada par de series
#'
#' @param control_data datos de las estaciones y sus coordenadas
#'
#' @return matriz con los meses que se solapan las estaciones entre si
#' @export
#'
overlap_station <- function(control_data){
  data_data <- control_data$data
  data_data[!is.na(data_data)] <- 1
  data_data[is.na(data_data)] <- 0

  overlap <- array(NA, dim = c(dim(data_data)[2], dim(data_data)[2]), dimnames = list(colnames(data_data), colnames(data_data)))

  stations <- colnames(data_data)
  station <- stations[1]
  for(station in stations){
    aux_overlap <- data_data[, station] + data_data
    aux_overlap <- aux_overlap >= 2
    overlap[station, ] <- apply(aux_overlap, c(2), sum)
  }
  return(overlap)
}

#' Rellenado mensual de las series
#' Usamos estaciones a menos de 200km con correlación por encima de 0.7
#' Para junio, julio y agosto, rellenamos con la más cercana
#' Utilizar el método que mejor correlaciona con la serie original
#'
#' @param control_data datos de las estaciones y sus coordenadas
#' @param min_correlation Correlación mínima para usar el dato en el relleno
#' @param max_dist distancia máxima entre las series a utilizar
#'
#' @return data y coor con los datos que no pasan el control eliminados
#' @export
#'
fill_series <- function(control_data, min_correlation, max_dist){

  return_control_data <- control_data

  ##  Para junio, julio y agosto, rellenamos con la más cercana
  data_near <- near_estations(data = control_data$data, coor = control_data$coor, max_dist = max_dist)

  all_series <- colnames(control_data$data)
  summer = grepl("Jun", rownames(return_control_data$data)) | grepl("Jul", rownames(return_control_data$data)) | grepl("Aug", rownames(return_control_data$data))

  i_serie <- all_series[1]
  for(i_serie in all_series){
    i_order <- data_near[[i_serie]] #Orden de distancia con las otras series
    data_10 <- control_data$data[, i_order] # Datos en orden de las otras series
    data_refe_10 <- t(apply(data_10, c(1), select_data, n_reference_stations = 1))
    return_control_data$data[summer & is.na(return_control_data$data[, i_serie]), i_serie] <- data_refe_10[summer & is.na(return_control_data$data[, i_serie])]
  }

  ## Rellenamos con estaciones a menos de 200km con correlación por encima de 0.7
  # Las series a utilizar se tienen que solapar 10 años
  overlap <- overlap_station(control_data)

  # Rellenamos los meses que no son de verano
  data_cor_list <- near_correlations(data = control_data$data, coor = control_data$coor, max_dist = max_dist)
  months <- names(data_cor_list)
  months <- months[! months %in% c("Jun", "Jul", "Aug")]
  i_month <- months[2]
  for(i_month in months){
    data_cor <- data_cor_list[[i_month]]
    date_month <- grepl(i_month, rownames(control_data$data))
    data_month <- control_data$data[date_month, ]

    # No rellenamos con datos que no cumplan la correlación mínima
    data_cor[data_cor < min_correlation] <- NA

    # No rellenamos con datos de estaciones sin un solapamiento mínimo
    data_cor[overlap < min_overlap * 12] <- NA

    #Ordenar según correlaciones
    if(length(data_near) > 0){
      i_series <- all_series[1]
      for(i_series in all_series){
        i_order <- all_series[data_near[[i_series]]] #Estaciones ordenadas según su cercanía a la estación
        data_cor_series <- colnames(data_cor)[!is.na(data_cor[, i_series])]
        select_series <- i_order[i_order %in% data_cor_series] #Eliminamos estaciones que no cumplen alguna restricción
        other_series <- data_month[, select_series]
        if(is.null(dim(other_series))){
          other_series <- array(other_series, c(length(other_series), 1))
        }

        data_refe <- fill_one_series(series = data_month[, i_series], other_series = other_series)

        return_control_data$data[date_month & is.na(return_control_data$data[, i_series]), i_series] <- data_refe[is.na(return_control_data$data[date_month, i_series])]
      }
    }
  }
  return(return_control_data)
}

#' Rellena la serie recibida utilizando las otras en el orden en el que están en other_series
#'
#' @param series serie de datos a completar 
#' @param other_series series de datos con las que copletar en el orden en el que se tienen que utilizar 
#'
#' @return serie de datosd rellena
#' @export
#'
fill_one_series <- function(series, other_series){

  return_series <- series

  # Por cernaniá, solo miramos correlación para eliminar, 5 años comunes

  i_other_series = 1
  while(sum(is.na(return_series))> 0 & i_other_series <= dim(other_series)[2]){

    use_series <- other_series[,i_other_series]

    # Datos NA que son 0 en la serie de relleno se ponen directamente a 0
    return_series[is.na(return_series) & !is.na(use_series) & use_series == 0] = 0

    use_series_no_0 <- use_series[(is.na(use_series) | use_series != 0) & (is.na(series) | series != 0)]
    series_no_0 <- series[(is.na(use_series) | use_series != 0) & (is.na(series) | series != 0)]

    use_series_comoon = use_series_no_0[!is.na(series_no_0) & !is.na(use_series_no_0)]
    series_comoon = series_no_0[!is.na(series_no_0) & !is.na(use_series_no_0)]

    # Calculamos sobre los datos en común de la serie a usar para rellenar
    l_mom <- lmom::samlmu(use_series_comoon, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)

    # Seleccionamos la mejor de las 10 funciones
    p.value = array(NA, dim = length(pel_functions), dimnames = list(names(pel_functions)))
    best_funcion = names(pel_functions)[3]
    for(best_funcion in names(pel_functions)){
      shapiro_value <- tryCatch({
        par_exp <- pel_functions[[best_funcion]](l_mom)
        cdf_exp <- cdf_functions[[best_funcion]](use_series_comoon, para = par_exp)
        norm_exp <- stats::qnorm(cdf_exp)
        shapiro_value <- stats::shapiro.test(norm_exp)$p.value
      }, error = function(cond) {
        NA
      })
      p.value[best_funcion] <- shapiro_value
    }
    if(sum(!is.na(p.value)) > 0){
      best_funcion <- names(pel_functions)[which(p.value == max(p.value, na.rm = TRUE))[1]]

      # Utlizamos la mejor de las 10 funciones
      par_exp <- pel_functions[[best_funcion]](l_mom)

      # Calculamos sobre todos los datos de la serie a usar para rellenar
      cdf_exp <- cdf_functions[[best_funcion]](use_series_no_0, para = par_exp)

      generate_series <- as.vector(stats::quantile(series_no_0, cdf_exp, na.rm = TRUE))
      names(generate_series) <- names(series_no_0)
      return_series[names(return_series)[is.na(return_series)]] <- generate_series[names(return_series)[is.na(return_series)]]
    }
    i_other_series <- i_other_series + 1
  }

  return(return_series)
}

#' Guarda la salida en 5 fichero con los datos
#'  5 ficheros que indican si cada datos es original o rellenado (1 dato no alterado, 0 dato alterado)
#'  y 5 ficheros de coordenadas para las estaciones de cada fichero de datos, que son:
#' - 1870 a 2020 con más de 80 años originales
#' - 1900 a 2020 con más de 80 años originales
#' - 1930 a 2020 con más de 60 años originales
#' - 1950 a 2020 con más de 40 años originales
#' - 1990 a 2020 con más de 30 años originales
#'
#' @param data_ori datos originales leidos de los ficheros CSV
#' @param control_data datos de las estaciones y sus coordenadas
#' @param folder_name carpeta en la que se guardan los datos, si NA no se guardan
#'
#' @return data y coor con los datos que no pasan el control eliminados
#' @export
#'
save_data <- function(data_ori, control_data, folder_name = NA){
  data_return <- list()

  i_ini <- i_inis[1]
  for(i_ini in i_inis){
    i_year <- i_years[as.character(i_ini)]
    year_ini <- max(i_ini, as.numeric(substr(rownames(data_ori)[1], 8, 11)))
    year_end <- min(i_end, as.numeric(substr(rownames(data_ori)[length(rownames(data_ori))], 8, 11)))
    use_data_ini <- which(grepl(year_ini, rownames(data_ori)))
    use_data_end <- which(grepl(year_end, rownames(data_ori)))
    use_data <- use_data_ini[1]:use_data_end[length(use_data_end)]
    data <- data_ori[use_data, ]
    data <- data[, apply(data, c(2), sum_no_nas) > i_year * 12]
    coor_save <- control_data$coor[colnames(data), ]
    data_save <- control_data$data[rownames(control_data$data) %in% rownames(data), colnames(data)]

    data_return[[paste0("start_", i_ini)]] <- list(data_save = data_save, coor_save = coor_save)

    if(dim(data_save)[2] > 0 & !is.na(folder_name)){
      dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)

      save_csvs(i_ini, folder_name, data_save, coor_save)

      # Indicamos con un 1 que no se ha reconstruido el dato y con un 0 que si
      data_ori_save <- data_ori[rownames(data_save), colnames(data_save)]
      data_save_change <- data_save
      data_save_change[, ] <- 1
      data_save_change[data_ori_save != data_save] <- 0
      data_save_change[is.na(data_ori_save) & !is.na(data_save)] <- 0

      data_save_year_month <- as.data.frame(data_save_change)
      data_save_year_month["year"] = substr(rownames(data_save_year_month), 8, 11)
      data_save_year_month["month"] = as.numeric(months(chron::chron(rownames(data_save_year_month), format=c(dates = "d/m/yy", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s"))))
      data_save_year_month = data_save_year_month[c("year", "month", colnames(data_save_year_month)[1:(dim(data_save_year_month)[2]-2)])]

      utils::write.table(data_save_year_month, file.path(folder_name, paste0("data_change_", i_ini, ".csv")), row.names = FALSE, sep = ";")
    }
  }
  return(data_return)
}

#' Guardamos los datos en CSVs
#'
#' @param i_ini identificador de los ficheros
#' @param folder_name carpeta en la que guardar el fichero
#' @param data_save datos de las estaciones a guardar
#' @param coor_save datos de coordenadas a guardar
#'
#' @return None
#' @export
#'
save_csvs <- function(i_ini, folder_name, data_save, coor_save){
  if(length(data_save) > 0){
    dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)
    data_save_year_month <- as.data.frame(data_save)
    data_save_year_month["year"] = substr(rownames(data_save_year_month), 8, 11)
    data_save_year_month["month"] = as.numeric(months(chron::chron(rownames(data_save_year_month), format=c(dates = "d/m/yy", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s"))))
    data_save_year_month = data_save_year_month[c("year", "month", colnames(data_save_year_month)[1:(dim(data_save_year_month)[2]-2)])]

    utils::write.table(data_save_year_month, file.path(folder_name, paste0("data_", i_ini, ".csv")), row.names = FALSE, sep = ";")

    coor_save_csv <- as.data.frame(coor_save)
    coor_save_csv[, "coor"] <- rownames(coor_save_csv)
    coor_save_csv <- coor_save_csv[, c("coor", "lat", "lon")]
    utils::write.table(coor_save_csv, file.path(folder_name, paste0("coor_", i_ini, ".csv")), sep = ";", row.names = FALSE)
  }
}

#' Guardamos las variables globales utilizadas
#'
#' @param folder_name carpeta en la que guardar el fichero
#' @param max_dist máxima distancia entre 2 estaciones para usar una para evaluar o completar la otra
#'
#' @return None
#' @export
#'
save_parameters <- function(folder_name, max_dist){
  text <- paste("max_dist <-",  max_dist, "
  min_years <-", min_years, "
  min_overlap <-", min_overlap, "
  n_reference_stations <-", n_reference_stations, "
  min_correlation <-", min_correlation, "
  max_diff_anomaly <-", max_diff_anomaly, "
  max_diff_anomaly_0 <-", max_diff_anomaly_0, "
  percentage_filled_data <-", percentage_filled_data, "
  min_second_correlation <-", min_second_correlation, "
  ")
  writeLines(text, file.path(folder_name, "parameters.txt"))
}

#' Leemos los datos desde los CSVs con el formato acordado
#'
#' @param file_data ruta del fichero de datos
#' @param file_coor ruta del fichero de coordenadas
#'
#' @return None
#' @export
#'
read_data <- function(file_data, file_coor){
  # Leer los datos
  if(file.exists(file_data) & file.exists(file_coor)){
    data_ori <- utils::read.table(file_data, sep = ";", header = TRUE)
    coor <- utils::read.table(file_coor, sep = ";", header = TRUE)
    rownames(coor) <- coor[, "coor"]
    coor[, "coor"] <- NULL

    # Las fechas identifican las filas
    rownames(data_ori) <- as.character(chron::chron(paste0("01", "/", data_ori[, "month"], "/", data_ori[, "year"]), format = c(dates = "d/m/y", times = "h:m:s"), out.format = time_format))
    data_ori["month"] <- NULL
    data_ori["year"] <- NULL

    # Revisión de fechas
    dates <- seq(chron::chron(rownames(data_ori)[1], format = time_format, out.format = time_format), chron::chron(rownames(data_ori)[dim(data_ori)[1]], format = time_format, out.format = time_format), by = "month")
    data <- data_ori[as.character(dates), ]

    if(is.null(dim(data))){
      data <- array(data, c(length(data), 1), list(as.character(dates), c(colnames(data_ori))))
    }

    # Revisión de que tenemos las coordenadas
    no_coor <- which(!colnames(data) %in% rownames(coor))
    if(length(no_coor) > 0){
      warning(paste("Eliminamos las estaciones", paste(colnames(data)[no_coor], sep = " ", collapse = " "), "por no tener disponibles las coordenadas.", sep = " ", collapse = " "))
      data <- data[, colnames(data) %in% rownames(coor)]
    }

    # Seleccionamos coordenadas solo de las estaciones con datos
    coor <- coor[colnames(data), ]

    # Eliminamos datos de fechas que no nos interesan
    start_date <- "01/01/1870"
    if(sum(dates == chron::chron(start_date)) > 0){
      dates <- dates[which(dates == chron::chron(start_date)):length(dates)]
      data <- data[as.character(dates), ]
    }
    return(list(data_ori = data_ori, data = data, coor = coor))
  }else{
    return(NA)
  }
}

# https://drive.google.com/file/d/1IvsjiOWdPUFzvM-DqDNYiLKRzAj7elwG/view
# Function snht - SNHT for (one) shift without trend
# x = test series Q with 1. column: year and 2. column: value
# critical level = c(90,92,94,95,97.5,99)
snht <- function(x, clevel)
{
  if(is.na(match(clevel,c(90, 92, 94, 95, 97.5, 99)))){
     stop("clevel must be 90, 92, 94, 95, 97.5 or 99")
  }
  n <- NROW(x)
  x1 <- x[, 1]
  y1 <- x[, 2]
  Z <- scale(y1)[, 1]
  Tv <- c()
  mean.1.2 <- array(data = NA, dim = c(0, 2))

  # critical levels (simulated by Khaliq and Ouarda 2007)
  cval <- matrix(data = c(
    10, 4.964, 5.197, 5.473, 5.637, 6.188, 6.769, 
    12, 5.288, 5.554, 5.876, 6.068, 6.729, 7.459, 
    14, 5.54, 5.831, 6.187, 6.402, 7.152, 8.001, 
    16, 5.749, 6.059, 6.441, 6.674, 7.492, 8.44, 
    18, 5.922, 6.248, 6.652, 6.899, 7.775, 8.807, 
    20, 6.07, 6.41, 6.83, 7.089, 8.013, 9.113, 
    22, 6.2, 6.551, 6.988, 7.257, 8.22, 9.38, 
    24, 6.315, 6.675, 7.123, 7.4, 8.4, 9.609, 
    26, 6.417, 6.785, 7.246, 7.529, 8.558, 9.812, 
    28, 6.509, 6.884, 7.353, 7.643, 8.697, 9.993, 
    30, 6.592, 6.973, 7.451, 7.747, 8.825, 10.153, 
    32, 6.669, 7.056, 7.541, 7.841, 8.941, 10.3, 
    34, 6.741, 7.132, 7.625, 7.93, 9.05, 10.434, 
    36, 6.803, 7.201, 7.699, 8.009, 9.143, 10.552, 
    38, 6.864, 7.263, 7.768, 8.081, 9.23, 10.663, 
    40, 6.921, 7.324, 7.835, 8.151, 9.317, 10.771, 
    42, 6.972, 7.38, 7.894, 8.214, 9.39, 10.865, 
    44, 7.022, 7.433, 7.951, 8.273, 9.463, 10.957, 
    46, 7.071, 7.484, 8.007, 8.331, 9.53, 11.04, 
    48, 7.112, 7.529, 8.054, 8.382, 9.592, 11.116, 
    50, 7.154, 7.573, 8.103, 8.432, 9.653, 11.193, 
    52, 7.194, 7.616, 8.149, 8.48, 9.711, 11.259, 
    54, 7.229, 7.654, 8.19, 8.524, 9.76, 11.324, 
    56, 7.264, 7.69, 8.23, 8.566, 9.81, 11.382, 
    58, 7.299, 7.727, 8.268, 8.606, 9.859, 11.446, 
    60, 7.333, 7.764, 8.308, 8.647, 9.906, 11.498, 
    62, 7.363, 7.796, 8.343, 8.683, 9.947, 11.548, 
    64, 7.392, 7.827, 8.375, 8.717, 9.985, 11.599, 
    66, 7.421, 7.857, 8.408, 8.752, 10.026, 11.648, 
    68, 7.449, 7.886, 8.439, 8.784, 10.067, 11.692, 
    70, 7.475, 7.913, 8.467, 8.814, 10.099, 11.737, 
    72, 7.499, 7.938, 8.496, 8.844, 10.134, 11.776, 
    74, 7.525, 7.965, 8.523, 8.873, 10.171, 11.822, 
    76, 7.547, 7.989, 8.548, 8.898, 10.2, 11.858, 
    78, 7.57, 8.013, 8.575, 8.926, 10.23, 11.895, 
    80, 7.591, 8.035, 8.599, 8.951, 10.259, 11.928, 
    82, 7.613, 8.059, 8.623, 8.976, 10.29, 11.966, 
    84, 7.634, 8.079, 8.647, 9.001, 10.315, 11.995, 
    86, 7.655, 8.102, 8.67, 9.026, 10.347, 12.033, 
    88, 7.673, 8.121, 8.691, 9.047, 10.37, 12.059, 
    90, 7.692, 8.14, 8.71, 9.067, 10.394, 12.089, 
    92, 7.711, 8.16, 8.732, 9.09, 10.417, 12.12, 
    94, 7.73, 8.181, 8.752, 9.11, 10.447, 12.153, 
    96, 7.745, 8.196, 8.77, 9.127, 10.465, 12.175, 
    98, 7.762, 8.214, 8.788, 9.147, 10.484, 12.196, 
    100, 7.778, 8.231, 8.807, 9.167, 10.507, 12.228, 
    105, 7.819, 8.273, 8.851, 9.214, 10.562, 12.291, 
    110, 7.856, 8.312, 8.892, 9.255, 10.608, 12.343, 
    115, 7.891, 8.35, 8.931, 9.296, 10.656, 12.401, 
    120, 7.921, 8.38, 8.963, 9.33, 10.694, 12.446, 
    125, 7.952, 8.413, 8.999, 9.365, 10.735, 12.488, 
    130, 7.983, 8.446, 9.032, 9.4, 10.772, 12.538, 
    135, 8.01, 8.474, 9.063, 9.431, 10.808, 12.579, 
    140, 8.038, 8.501, 9.092, 9.462, 10.845, 12.621, 
    145, 8.063, 8.529, 9.12, 9.49, 10.877, 12.66, 
    150, 8.086, 8.554, 9.147, 9.519, 10.906, 12.694, 
    155, 8.111, 8.578, 9.172, 9.543, 10.933, 12.725, 
    160, 8.133, 8.601, 9.195, 9.569, 10.966, 12.759, 
    165, 8.155, 8.625, 9.222, 9.596, 10.992, 12.793, 
    170, 8.174, 8.643, 9.241, 9.615, 11.016, 12.82, 
    175, 8.195, 8.666, 9.265, 9.641, 11.046, 12.851, 
    180, 8.214, 8.685, 9.283, 9.658, 11.062, 12.872, 
    185, 8.233, 8.706, 9.307, 9.683, 11.089, 12.904, 
    190, 8.252, 8.725, 9.325, 9.701, 11.11, 12.93, 
    195, 8.268, 8.741, 9.343, 9.72, 11.132, 12.956, 
    200, 8.286, 8.761, 9.364, 9.741, 11.156, 12.982, 
    225, 8.361, 8.838, 9.446, 9.826, 11.247, 13.083, 
    250, 8.429, 8.908, 9.516, 9.898, 11.329, 13.175, 
    275, 8.489, 8.97, 9.581, 9.966, 11.399, 13.248, 
    300, 8.54, 9.022, 9.635, 10.02, 11.46, 13.326, 
    325, 8.587, 9.07, 9.685, 10.071, 11.517, 13.389, 
    350, 8.633, 9.117, 9.732, 10.118, 11.565, 13.44, 
    375, 8.67, 9.157, 9.775, 10.161, 11.613, 13.494, 
    400, 8.706, 9.193, 9.814, 10.202, 11.654, 13.542, 
    425, 8.738, 9.224, 9.844, 10.234, 11.692, 13.58, 
    450, 8.771, 9.26, 9.882, 10.272, 11.73, 13.623, 
    475, 8.798, 9.288, 9.912, 10.302, 11.761, 13.655, 
    500, 8.828, 9.317, 9.939, 10.33, 11.795, 13.69, 
    525, 8.854, 9.344, 9.967, 10.36, 11.827, 13.73, 
    550, 8.878, 9.369, 9.995, 10.386, 11.854, 13.751, 
    575, 8.901, 9.391, 10.016, 10.408, 11.878, 13.782, 
    600, 8.923, 9.414, 10.04, 10.431, 11.904, 13.813, 
    650, 8.963, 9.455, 10.083, 10.476, 11.949, 13.856, 
    700, 9.001, 9.493, 10.119, 10.511, 11.986, 13.904, 
    750, 9.033, 9.524, 10.152, 10.547, 12.026, 13.947, 
    800, 9.063, 9.557, 10.187, 10.58, 12.059, 13.975, 
    850, 9.093, 9.587, 10.216, 10.612, 12.096, 14.023, 
    900, 9.119, 9.614, 10.244, 10.64, 12.12, 14.041, 
    950, 9.143, 9.638, 10.269, 10.665, 12.149, 14.07, 
    1000, 9.168, 9.664, 10.295, 10.692, 12.176, 14.105, 
    1100, 9.211, 9.708, 10.339, 10.736, 12.22, 14.15, 
    1200, 9.246, 9.745, 10.377, 10.775, 12.263, 14.197, 
    1300, 9.283, 9.781, 10.415, 10.812, 12.304, 14.235, 
    1400, 9.313, 9.812, 10.446, 10.845, 12.34, 14.271, 
    1500, 9.347, 9.846, 10.481, 10.88, 12.374, 14.312, 
    1600, 9.372, 9.871, 10.506, 10.904, 12.396, 14.339, 
    2000, 9.464, 9.965, 10.603, 11.002, 12.5, 14.443, 
    2500, 9.551, 10.052, 10.69, 11.089, 12.591, 14.54, 
    3000, 9.618, 10.121, 10.76, 11.161, 12.664, 14.619, 
    3500, 9.675, 10.178, 10.818, 11.219, 12.727, 14.683, 
    4000, 9.727, 10.229, 10.869, 11.271, 12.779, 14.734, 
    4500, 9.766, 10.269, 10.911, 11.313, 12.82, 14.777, 
    5000, 9.803, 10.307, 10.948, 11.349, 12.859, 14.817, 
    7500, 9.938, 10.442, 11.085, 11.487, 12.997, 14.959, 
    10000, 10.031, 10.537, 11.18, 11.584, 13.095, 15.063, 
    15000, 10.152, 10.658, 11.302, 11.707, 13.221, 15.186, 
    20000, 10.236, 10.743, 11.388, 11.791, 13.305, 15.271, 
    50000, 10.48, 10.988, 11.634, 12.039, 13.556, 15.523), nrow = 108, ncol = 7, byrow = TRUE,
    dimnames = list(NULL, c('n', 'clevel90', 'clevel92', 'clevel94', 'clevel95', 'clevel97.5', 'clevel99')))

  for (v in 1:(n-1)) {
    z1 <- mean(Z[1:v], na.rm=T)
    z2 <- mean(Z[(v + 1):n], na.rm=T)
    mean.1.2 <- rbind(mean.1.2, c(z1, z2))
    Tvi <- (v * (z1^2)) + ((n - v) * (z2^2))
    Tv <- rbind(Tv, Tvi)
  }
  # Test Statistic (omit tail ends +-5 years)
  T0 <- max(Tv[5:(length(Tv) - 5)])
  T0x <- which(Tv == T0)
  T0xa <- as.numeric(x1[which(Tv == T0)])
  Tc <- as.numeric(cval[max(which(cval[, 'n'] <= n)), paste('clevel', clevel, sep = "")])
  if (T0 > Tc) {
    significant <- 1
  } else {
    significant <- 0
  }
  mean.1 <- mean.1.2[T0x, 1] * stats::sd(y1, na.rm=T) + mean(y1, na.rm=T)
  mean.2 <- mean.1.2[T0x, 2] * stats::sd(y1, na.rm=T) + mean(y1, na.rm=T)

  res <- list(x = x, Tv = Tv, Tc = Tc, T0 = T0, T0x = T0x, T0xa = T0xa,
              mean.1 = mean.1, mean.2 = mean.2,
              clevel = clevel, sig = significant)
  return(res)
}

#' Devuelve el pval calculado por mkTrend o el pval0 si el pval era NA
#'
#' @param data matriz de datos
#'
#' @return pval
#' @export
#'
calc_mkTrend_pval <- function(data) {
  mkTrend_data <- mkTrend(data)
  mkTrend_pval <- mkTrend_data$`Corrected p.value`
  mkTrend_pval0 <- mkTrend_data$p.value
  mkTrend_pval[is.na(mkTrend_pval) | is.null(mkTrend_pval)] <- mkTrend_pval0[is.na(mkTrend_pval) | is.null(mkTrend_pval)]

  if(sum(data != 0) == 0){
    mkTrend_pval[1] = 1
  }

  return(mkTrend_pval)
}

#' Suma los datos de cada año, para devolver un solo dato anual
#'
#' @param data matriz de datos
#'
#' @return un dato por año
#' @export
#'
calc_data_year <- function(data) {
  year <- as.numeric(substr(rownames(data), 8, 11))
  data_year <- as.data.frame(data) 
  data_year <- stats::aggregate(data_year, by = list("year" = year), FUN = sum)
  rownames(data_year) <- paste0("year_", as.character(unique(year)))
  data_year$year <- NULL
  return(as.matrix(data_year))
}

#' Devuelve el slope z por meses, años y estaciones
#'
#' @param data datos de las estaciones
#' @param calc_function función a utilizar
#'
#' @return lista de resultados
#' @export
#'
calc_data_year_month_station <- function(data, calc_function) {
  months <- substr(rownames(data), 4, 6)

  # Tendencia mensual, estacional y anual,función calc_function
  sens_slope <- list()
  sens_slope[["year"]] <- apply(calc_data_year(data), c(2), calc_function)
  sens_slope[["summer"]] <- apply(calc_data_year(data[months %in% c("Jun", "Jul", "Aug"), ]), c(2), calc_function)
  sens_slope[["autumn"]] <- apply(calc_data_year(data[months %in% c("Sep", "Oct", "Nov"), ]), c(2), calc_function)
  sens_slope[["winter"]]  <- apply(calc_data_year(data[months %in% c("Dec", "Jan", "Feb"), ]), c(2), calc_function)
  sens_slope[["spring"]] <- apply(calc_data_year(data[months %in% c("Mar", "Apr", "May"), ]), c(2), calc_function)
  month <- months[8]
  for(month in unique(months)){
    sens_slope[[month]] <- apply(calc_data_year(data[months == month, ]), c(2), calc_function)
  }
  return(sens_slope)

  # i <- 12
  # trend <- list()
  # trend[["year"]] <- data[, i]
  # trend[["summer"]] <- data[months %in% c("Jun", "Jul", "Aug"), i]
  # trend[["autumn"]] <- data[months %in% c("Sep", "Oct", "Nov"), i]
  # trend[["winter"]]  <- data[months %in% c("Dec", "Jan", "Feb"), i]
  # trend[["spring"]] <- data[months %in% c("Mar", "Apr", "May"), i]
  # month <- months[8]
  # for(month in unique(months)){
  #   trend[[month]] <- data[months == month, i] 
  # }
  # saveRDS(trend, "trend12.rds")
}

#' esta función calcula las tendencias. Lo que nos interesa es slp (es la función de cálculo del Sen slope), pval (a veces no da resultado por temas de iteración) entonces coger pval0.
#'
#' @param x x
#' @param ci ci
#'
#' @return list
#' @export
#'
mkTrend <- function(x, ci = .95) {
  x <- x + 1 # Añadido para evitar algunos problemas con los 0s

  z <- NULL
  z0 <- NULL
  pval <- NULL
  pval0 <- NULL
  S <- 0
  Tau <- NULL
  essf <- NULL
  if (!is.vector(x)) {
    stop("Input data must be a vector")
  }
  if (any(!is.finite(x))) {
    x <- x[- c(which(!is.finite(x)))]
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }
  n <- length(x)
  for (i in 1:(n - 1)) {
   for (j in (i + 1):n) {
    S <- S + sign(x[j] - x[i])
   }
  }
  ro <- stats::acf(rank(stats::lm(x ~ I(1:n))$resid), lag.max = (n - 1), plot = FALSE)$acf[-1]
  sig <- stats::qnorm((1 + ci) / 2) / sqrt(n)
  rof <- rep(NA, length(ro)) 
  for (i in 1:(length(ro))) {
   if(ro[i] > sig || ro[i] < -sig) {
    rof[i] <- ro[i]
   } else {
    rof[i] <- 0
   }
  }
  cte <- 2 / (n * (n - 1) * (n - 2))
  ess <- 0
  for (i in 1:(n-1)) {
    ess <- ess + (n - i) * (n - i - 1) * (n - i - 2) * rof[i]
  }
  essf <- 1 + ess * cte
  var.S <- n * (n - 1) * (2 * n + 5) * (1 / 18)
  if(length(unique(x)) < n) {
   aux <- unique(x) 
   for (i in 1:length(aux)) {
     tie <- length(which(x == aux[i]))
     if (tie > 1) {
      var.S <- var.S - tie * (tie - 1) * (2 * tie + 5) * (1 / 18)
     }
   }
  }
  VS <- var.S * essf
  if (S == 0) {
   z <- 0
   z0 <- 0
  }
  if (S > 0) {
   z <- (S - 1) / sqrt(VS)
   z0 <- (S - 1)/sqrt(var.S)
  } else {
   z <- (S + 1) / sqrt(VS)
   z0 <- (S + 1) / sqrt(var.S)
  }
  pval <- 2 * stats::pnorm(- abs(z))
  pval0 <- 2 * stats::pnorm(- abs(z0))
  Tau <- S / (.5 * n * (n - 1))
  V <- rep(NA, times = (n^2 - n) / 2)
  k <- 0
  for (i in 2:n) {
   for (j in 1:(n - 1)) {
     k <- k + 1
     V[k] <- (x[i] - x[j]) / (i - j)
   }
  }
  slp <- stats::median(stats::na.omit(V))
  return(list("Z" = z0, "p.value" = pval0, "Zc" = z, "Corrected p.value" = pval, "tau" = Tau, "N/N*s" = essf, "Sen's Slope" = slp))
}

#' regresión lineal de los datos contra los años
#'
#' @param data index
#'
#' @return lm
#' @export
#'
calc_mkTrend_slp <- function(data){
  years <- as.numeric(substr(names(data), 8, 11))
  # mkTrend_slp <- stats::lm(data ~ years)$coefficients[[1]]
  if(sum(data) != 0){
    # mkTrend_slp <- mkTrend(data)$"Sen's Slope"
    mkTrend_slp <- RobustLinearReg::theil_sen_regression(data ~ years)$coefficients[2]
  }else{
    mkTrend_slp <- 0  
  }
  return(mkTrend_slp)
}

#############################Funciones de Sergio

#' esta función calcula la tendencia. Hay que definirle un objeto de años (years) con el año correspondiente a cada caso.
#'
#' @param index index
#' @param threshold threshold
#'
#' @return output
#' @export
#'
dry_spell_trend <- function(index, threshold) {
  a.1 <- NA
  a.2 <- NA
  a.3 <- NA
  a.4 <- NA

  if(sum(!is.na(index)) > 0){
    years <- as.numeric(substr(names(index), 8, 11))
    n <- length(index)
    below <- which(index < threshold)
    above <- which(index >= threshold)
    dry <- index
    dry[below] <- 1
    dry[above] <- NA
    index2 <- - dry * index
    duration1 <- rep(NA,n)
    magnitud1 <- rep(NA,n)
    for (i in 2:n){
      if(!is.na(dry[i])){
        if(is.na(dry[i - 1]) & dry[i] == 1){
          duration1[i] <- dry[i]
        }else{
          if(dry[i - 1] == 1 & dry[i] == 1){
            duration1[i] <- duration1[i - 1] + dry[i]
          }
        }
      }
      if(is.na(index2[i])){
        magnitud1[i]
      }else{
        if(is.na(index2[i - 1]) & !is.na(index2[i])){
          magnitud1[i] <- index2[i]
        }else{
          if(!is.na(index2[i - 1]) & !is.na(index2[i])){
            magnitud1[i] <- magnitud1[i - 1] + index2[i]
          }else{
            magnitud1[i] <- NA
          }
        }
      }
    }
    duration2 <- rep(NA, n)
    magnitud2 <- rep(NA, n)
    for (i in 1:n){
      if(!is.na(duration1[i]) & is.na(duration1[i + 1])){
        duration2[i] <- duration1[i]
      }    
      if(!is.na(magnitud1[i]) & is.na(magnitud1[i + 1])){
        magnitud2[i] <- magnitud1[i]
      }      
    }
    cases <- which(!is.na(duration2))
    series_dry_ori <- cbind(years[cases], duration2[cases], magnitud2[cases])
    colnames(series_dry_ori) <- c("year", "duration", "magnitude")
    series_dry <- array(0, dim = c(length(unique(years)), 3), dimnames = list(as.character(unique(years)), colnames(series_dry_ori)))
    series_dry[, "year"] <- as.numeric(rownames(series_dry))
    series_dry_ori <- as.data.frame(series_dry_ori)

    if(dim(series_dry_ori)[1] > 0){
      series_dry_aux <- as.matrix(stats::aggregate(series_dry_ori, by = list("year" = series_dry_ori$year), FUN = sum))
      series_dry[as.character(series_dry_aux[, c("year")]), c("duration", "magnitude")] <- as.numeric(series_dry_aux[, c("duration", "magnitude")])
    }    
    if(dim(series_dry_ori)[1] > 2){
      # a.1 <- as.vector(Kendall::Kendall(series_dry[,1], series_dry[,2])$sl)
      a.1 <- calc_mkTrend_pval(series_dry[, 2])
      if(sum(is.infinite(series_dry)) <= 0){
        # a.2 <- as.vector(Kendall::Kendall(series_dry[,1], series_dry[,3])$sl)
        a.2 <- calc_mkTrend_pval(series_dry[, 3])
      }
    } 
    if(dim(series_dry_ori)[1] > 0){
      # lm.1 <- stats::lm(series_dry[,2] ~ series_dry[,1])
      # kk.1 <- as.data.frame(stats::coef(lm.1))[2,]
      # kk.1 <- mkTrend(series_dry[, 2])$"Sen's Slope"
      # a.3 <- kk.1 * 10
      theil_sen_regression_data2 <- as.vector(series_dry_ori[, 2])
      theil_sen_regression_data1 <- as.vector(series_dry_ori[, 1])
      a.3 <- RobustLinearReg::theil_sen_regression(theil_sen_regression_data2 ~ theil_sen_regression_data1)$coefficients[2]
      if(sum(is.infinite(series_dry)) <= 0){
        # lm.2 <- stats::lm(series_dry[,3] ~ series_dry[,1])
        # kk.2 <- as.data.frame(stats::coef(lm.2))[2,]
        # kk.2 <- mkTrend(series_dry[, 3])$"Sen's Slope"
        # a.4 <- kk.2 * 10
        theil_sen_regression_data3 <- as.vector(series_dry_ori[, 3])
        theil_sen_regression_data1 <- as.vector(series_dry_ori[, 1])
        a.4 <- RobustLinearReg::theil_sen_regression(theil_sen_regression_data3 ~ theil_sen_regression_data1)$coefficients[2]
      }
    }
  }
  output <- cbind(a.1, a.2, a.3, a.4)
  colnames(output) <- c("p.dur", "p.mag", "trend.dur", "trend.mag")
  return(output)
}

#' esto te calcula unas tendencias móviles de una serie, en este caso que empieza en 1851 y termina en 2018,. Habría que hacerlo para cada base de datos y estación.
#'
#' @param datos datos
#'
#' @return list
#' @export
#'
mobile_trends <- function(datos) {
  anos <- unique(as.numeric(substr(names(datos), 8, 11)))

  matriz_p <- NA
  matriz_s <- NA
  matriz_per <- NA
  if(length(anos) > 30){
    matriz_p <- matrix(NA, length(anos) - 30 + 1, length(anos) - 30 + 1)
    matriz_s <- matrix(NA, length(anos) - 30 + 1, length(anos) - 30 + 1)
    matriz_per <- matrix(NA, length(anos) - 30 + 1, length(anos) - 30 + 1)

    bb <- 1
    for (bb in c(1:dim(matriz_p)[1])){
      aa <- length(anos) - bb + 1
      rr <- aa - 30 + 1
      for (i in 1:rr){
        if(sum(datos[i:aa]) != 0){
          mkTrend_data <- mkTrend(datos[i:aa])
          matriz_p[bb, i] <- mkTrend_data$p.value
          # matriz_s[bb, i] <- as.data.frame(stats::lm(datos[i:aa] ~ anos[i:aa])$coef)[2, ]
          theil_sen_regression_datos <- datos[i:aa]
          theil_sen_regression_anos <- anos[i:aa]
          matriz_s[bb, i] <- RobustLinearReg::theil_sen_regression(theil_sen_regression_datos ~ theil_sen_regression_anos)$coefficients[2]
          # matriz_s[bb, i] <- mkTrend_data$"Sen's Slope"
          matriz_per[bb, i] <- calc_percentage(datos[i:aa])
        }else{
          matriz_p[bb, i] <- 1
          matriz_s[bb, i] <- 0
          matriz_per[bb, i] <- 0
        }
      }
    }
  }
  return(list(matriz_p = matriz_p, matriz_s = matriz_s, matriz_per = matriz_per))
}

#' Diferencia en porcentaje
#'
#' @param datos datos
#'
#' @return percentage
#' @export
#'
calc_percentage <- function(datos) {
  years <- as.numeric(substr(names(datos), 8, 11))
  intercept_t <- RobustLinearReg::theil_sen_regression(datos ~ years)$coefficients
  final <- length(datos) * intercept_t[2] + intercept_t[1]
  inicial <- 1 * intercept_t[2] + intercept_t[1]
  percentage <- 100 * final / inicial - 100
  if(is.na(percentage)){
    percentage <- 0
  }
  return(percentage)
}
