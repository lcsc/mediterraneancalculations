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
#' @importFrom sf st_sfc st_multipoint st_transform st_coordinates
#' @importFrom SpatialTools dist1
#' @importFrom lmom samlmu pelexp cdfexp pelgam cdfgam pelgev cdfgev pelglo cdfglo pelgpa cdfgpa pelgno cdfgno pelln3 cdfln3 pelnor cdfnor pelpe3 cdfpe3 pelwei cdfwei
#' @importFrom SPEI spi
#' @importFrom hydroGOF d mae pbias rmse
#' @importFrom snowfall sfInit sfExport sfLapply sfStop sfExportAll sfLibrary
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

### Variables globales
crs84 <- 4326 #https://spatialreference.org/ref/epsg/wgs-84/ "epsg:4326"
crs84m <- 4978  #https://epsg.io/4978 "EPSG:4978"

min_years <- 10 # Años mínimos de datos para que no se elimine una estación
min_overlap <- 10 # Años mínimos de datos de solape entre 2 series para rellenar una con la otra
n_reference_stations <- 10 # Número de datos usados para sacar la serie de referencia en el control de calidad
min_correlation <- 0.7 # Correlación mínima para usar el dato en el relleno
max_diff_anomaly <- 0.8 # Máxima diferencia de anomalías para mantener dato en el control
max_diff_anomaly_0 <- 0.8 # Máxima diferencia de anomalías para mantener dato en el control, si el dato es 0
percentage_filled_data <- 95 # Número de datos que tiene una estación para entrar en el segundo relleno
min_second_correlation <- 0.5 # Correlación mínima entre estaciones para usarse en el segundo relleno

####################

# library(rworldmap, include.only = "getMap")
# world = rworldmap::getMap(resolution = "high")

time_format <- c(dates = "d/m/yy", times = "h:m:s")

pel_functions <- list(exp = lmom::pelexp, gam = lmom::pelgam, gev = lmom::pelgev, glo = lmom::pelglo, gpa = lmom::pelgpa, gno = lmom::pelgno, ln3 = lmom::pelln3, nor = lmom::pelnor, pe3 = lmom::pelpe3, wei = lmom::pelwei)
cdf_functions <- list(exp = lmom::cdfexp, gam = lmom::cdfgam, gev = lmom::cdfgev, glo = lmom::cdfglo, gpa = lmom::cdfgpa, gno = lmom::cdfgno, ln3 = lmom::cdfln3, nor = lmom::cdfnor, pe3 = lmom::cdfpe3, wei = lmom::cdfwei)

alexanderson_folder <- "results"

i_inis <- c(1871, 1916, 1931, 1951, 1961, 1981)
i_end <- 2020

# snht critical levels (simulated by Khaliq and Ouarda 2007)
snht_cval <- matrix(data = c(
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

####################

#' Lee los años de cadenas de texto que terminan con los años
#'
#' @param txt texto o vector de textos
#'
#' @return list
#' @export
#'
read_years <- function(txt){
  return(as.numeric(substr(txt, nchar(txt) - 4 + 1, nchar(txt))))
}

#' Ordena los datos y devuelve una lista con el orden
#'
#' @param data datos
#'
#' @return list
#' @export
#'
order_data <- function(data){
  data_order <- apply(data, c(1), order, na.last = NA)
  if("matrix" %in% class(data_order)){
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

  l_mom <- suppressWarnings(lmom::samlmu(data, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0))

  # Seleccionamos la mejor de las 10 funciones
  p.value = array(NA, dim = length(pel_functions), dimnames = list(names(pel_functions)))
  best_funcion = names(pel_functions)[3]
  for(best_funcion in names(pel_functions)){
    shapiro_value <- tryCatch({
      par_exp <- suppressWarnings(pel_functions[[best_funcion]](l_mom))
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
  points_data <- sf::st_sfc(sf::st_multipoint(as.matrix(coor[, c("lon", "lat"), drop = FALSE])), crs = crs84)
  # plot(world)
  # plot(points_data, add = TRUE, col = "red")
  points_data <- sf::st_transform(points_data, crs = crs84m)
  coordinates <- as.matrix(coor)
  coordinates[, ] <- sf::st_coordinates(points_data)[, c("X", "Y")]


  # Eliminar estaciones a más de 200 km para los cálculos del control de calidad
  dist <- SpatialTools::dist1(coordinates)
  colnames(dist) <- rownames(dist) <- rownames(coor)
  if(!is.na(max_dist)){
    dist[dist > max_dist * 1000] <- NA
  }
  dist[dist == 0] <- NA

  data_cor_list <- list()
  # Calcular correlaciones (por meses)
  while(length(data_cor_list) <= 11){
    data_month <- data[(c(1:dim(data)[1])-1)%%12 == length(data_cor_list), , drop = FALSE]
    name_month <- substr(rownames(data_month)[1], 4, 6)
    data_cor <- suppressWarnings(stats::cor(data_month, use = 'pairwise.complete.obs'))
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
  points_data <- sf::st_sfc(sf::st_multipoint(as.matrix(coor[, c("lon", "lat"), drop = FALSE])), crs = crs84)
  # plot(world)
  # plot(points_data, add = TRUE, col = "red")
  points_data <- sf::st_transform(points_data, crs = crs84m)
  coordinates <- as.matrix(coor)
  coordinates[, ] <- sf::st_coordinates(points_data)[, c("X", "Y"), drop = FALSE]

  # Eliminar estaciones a más de 200 km para los cálculos del control de calidad
  dist <- SpatialTools::dist1(coordinates)
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
#' Usando las 10 más cercanas a menos de 200 km, desechar si promedio de percentil se diferencia en más de 0.6 
#' o en más de 0.5 para datos 0
#'
#' @param data datos  
#' @param coor coordenadas
#' @param max_dist máxima distancia entre 2 estaciones para usar una para evaluar o completar la otra
#' @param max_diff_anomaly máxima diferencia de anomalías para mantener dato en el control 
#' @param max_diff_anomaly_0 máxima diferencia de anomalías para mantener dato en el control, si el dato es 0
#'
#' @return data y coor con los datos que no pasan el control eliminados
#' @export
#'
quality_control <- function(data, coor, max_dist, max_diff_anomaly, max_diff_anomaly_0){

  max_usable_neighbours <- 10
  min_usable_neighbours <- 4

  # Eliminar estaciones con menos de min_years años de datos
  min_n_data <- min_years * 12
  no_nas <- apply(data, c(2), sum_no_nas) > min_n_data
  data <- data[, no_nas, drop = FALSE]
  coor <- coor[colnames(data), , drop = FALSE]

  # Eliminar estaciones con menos de 1 año de datos para algún mes
  no_nas <- NA
  i <- 1
  for(i in c(0:11)){
    month_no_nas <- apply(data[c(1:dim(data)[1])%%12 == i, , drop = FALSE], c(2), sum_no_nas) > 1
    if(length(no_nas) > 1){
      no_nas <- no_nas & month_no_nas
    }else{
      no_nas <- month_no_nas
    }
  }
  data <- data[, no_nas, drop = FALSE]
  coor <- coor[colnames(data), , drop = FALSE]

  all_series <- colnames(data)
  data_correct <- as.matrix(data) 

  #Calcular anomalías de cada dato dentro de cada estación
  data_percent <- apply(data, c(2), apply_ecdf)

  ### Correlación por meses
  data_near <- near_estations(data, coor, max_dist = max_dist)

  if(length(data_near) > 0){
    months <- unique(substr(rownames(data), 4, 6))
    # Para cada mes
    i_month <- months[12]
    for(i_month in months){
      # Seleccionamos los datos de mes
      data_percent_month <- data_percent[grepl(i_month, rownames(data_percent)), ]

      #Ordenar según distancias
      data_order <- data_near
      if(length(data_order) > 0){
        # 10 series más cercanas
        i_serie <- all_series[1]
        for(i_serie in all_series){
          i_order <- data_order[[i_serie]][1:max_usable_neighbours] #Orden de cercanía con las otras series
          i_order <- i_order[!is.na(i_order)]
          if(length(i_order) >= min_usable_neighbours){
            data_ori <- data[rownames(data_percent_month), i_serie]
            data_10 <- data_percent_month[, i_order] # Datos en orden de las otras series
            data_refe_10 <- t(apply(data_10, c(1), select_data, n_reference_stations = n_reference_stations))
            data_refe <- apply(data_refe_10, c(1), base::mean, na.rm = TRUE)
            data_correct[rownames(data_percent_month)[!is.na(data_ori) & data_ori != 0 & !is.na(data_percent_month[, i_serie]) & !is.na(data_refe) & (abs(data_percent_month[, i_serie] - data_refe) > max_diff_anomaly)], i_serie] <- NA
            data_correct[rownames(data_percent_month)[!is.na(data_ori) & data_ori == 0 & !is.na(data_percent_month[, i_serie]) & !is.na(data_refe) & (abs(data_percent_month[, i_serie] - data_refe) > max_diff_anomaly_0)], i_serie] <- NA
          }
        }
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

#' Calcula el tiempo de solape existente entre cada par de series sin contar 0s
#'
#' @param control_data datos de las estaciones y sus coordenadas
#'
#' @return matriz con los meses que se solapan las estaciones entre si
#' @export
#'
overlap_station_no_0 <- function(control_data){
  data_data <- control_data$data

  data_data[data_data == 0] <- NA

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

#' En los países que no salgan series, vamos a permitir que hasta tres años de datos se rellenen con la media. Es decir, pongamos que si para un periodo concreto 1900-2020 no salen series pero saldrían porque hay un máximo de tres años de datos (es decir 36 meses), rellenamos esos datos con el promedio de los 15 datos más cercanos en el tiempo. Por ejemplo, si es 1900, pues con la media de 1900-1915, si es 1915, pues con la media de 1907 a 1922. Siempre y cuando esos cinco años no estén entre 2015 y 2020 o en los cinco primeros años de las series, que entonces tiramos la serie pues podría afectar a las tendencias.
#' Si las series son las de 1981-2020, lo mismo, pero dejamos solamente dos años de datos perdidos.
#'
#' @param data datos de las estaciones que se intentarán rellenar
#' @param fillable_years años rellenables con la media mensual de la propia estación
#'
#' @return None
#' @export
#'
fill_unfillable_station  <- function(data, fillable_years){
  fillable_data <- data$data[, apply(is.na(data$data), c(2), sum) > 0 & apply(is.na(data$data), c(2), sum) <= fillable_years & apply(is.na(data$data[1:60, , drop = FALSE]), c(2), sum) == 0 & apply(is.na(data$data[(dim(data$data)[1]-60+1):dim(data$data)[1], , drop = FALSE]), c(2), sum) == 0, drop = FALSE]

  if(!is.null(dim(fillable_data))){
    station <- colnames(fillable_data)[1]
    for(station in colnames(fillable_data)){
      data_padding <- fillable_data[, station]
      data_na <- which(is.na(fillable_data[, station]))
      data_no_na <- which(!is.na(fillable_data[, station]))

      idata_na <- data_na[1]
      for(idata_na in data_na){
        idata_no_na <- data_no_na[idata_na %% 12 == data_no_na %% 12]
        # rellenamos con la media de los 10 datos más cercanos
        idata_select_no_na <- idata_no_na[order(abs(idata_no_na - idata_na))][1:10]

        fillable_data[idata_na, station] <- base::mean(data_padding[idata_select_no_na])
      }
    }
    data$data[, colnames(fillable_data)] <- fillable_data
  }
  return(data)
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

  array_colnames <- c("d", "mae", "pbias", "rmse")

  return_control_data <- control_data

  if(dim(control_data$data)[2] > 1){
    all_series <- colnames(control_data$data)

    ## Rellenamos con estaciones a menos de 200km con correlación por encima de 0.7
    # Las series a utilizar se tienen que solapar 10 años
    overlap <- overlap_station(control_data)
    overlap_no_0 <- overlap_station_no_0(control_data)

    ##  Para algunos casos rellenamos con la más cercana
    data_near <- near_estations(data = control_data$data, coor = control_data$coor, max_dist = max_dist)

    # Rellenamos los meses que no son de verano
    data_cor_list <- near_correlations(data = control_data$data, coor = control_data$coor, max_dist = max_dist)
    months <- names(data_cor_list)
    i_month <- months[1]
    for(i_month in months){
      data_cor <- data_cor_list[[i_month]]
      date_month <- grepl(i_month, rownames(control_data$data))
      data_month <- control_data$data[date_month, ]

      test_serie <- array(NA, dim(data_month)[1])
      names(test_serie) <- rownames(data_month)

      overlap_month <- overlap_station(control_data = list(data = data_month, coor = control_data$coor))
      overlap_no_0_month <- overlap_station_no_0(control_data = list(data = data_month, coor = control_data$coor))

      # No rellenamos con datos que no cumplan la correlación mínima
      data_cor[data_cor < min_correlation] <- NA

      # No rellenamos con datos de estaciones sin un solapamiento mínimo para todo el periodo
      data_cor[overlap < min_overlap * 12] <- NA
      data_cor[overlap_no_0 < min_overlap * 12] <- NA

      # No rellenamos con datos de estaciones sin un solapamiento mínimo para el mes actual
      data_cor[overlap_month < min_overlap] <- NA
      data_cor[overlap_no_0_month < min_overlap] <- NA

      #Ordenar según correlaciones
      if(length(data_near) > 0){
        i_series <- all_series[1]
        for(i_series in all_series){
          i_order <- all_series[data_near[[i_series]]] #Estaciones ordenadas según su cercanía a la estación
          data_cor_series <- colnames(data_cor)[!is.na(data_cor[, i_series])]
          select_series <- i_order[i_order %in% data_cor_series] #Eliminamos estaciones que no cumplen alguna restricción
          other_series <- data_month[, select_series, drop = FALSE]

          if(dim(other_series)[2] != 0){
            data_refe <- fill_one_series(series = data_month[, i_series], other_series = other_series)
            data_refe_test <- fill_one_series(series = test_serie, other_series = other_series)
            return_control_data$data[date_month & is.na(return_control_data$data[, i_series]), i_series] <- data_refe[is.na(return_control_data$data[date_month, i_series])]           
          }

          zeros_left <- sum(is.na(data_month[, i_series])) > 0 && sum(data_month[, i_series] == 0, na.rm = TRUE)/sum(!is.na(data_month[, i_series]), na.rm = TRUE) > 0.5
          if(dim(other_series)[2] == 0 || zeros_left){
            i_order <- all_series[data_near[[i_series]]] #Orden de distancia con las otras series
            
            if(length(i_order) > 1){
              data_10 <- control_data$data[, i_order] # Datos en orden de las otras series
              data_refe_10 <- t(apply(data_10, c(1), select_data, n_reference_stations = 1))

              # Cuando se sustituye por la más cercana, considerar la ratio entre la estación candidata y la que utilizamos para sustituir
              mean_10 <- mean(return_control_data$data[, i_series][!is.na(return_control_data$data[, i_series]) & !is.na(data_refe_10)])
              mean_refe_10 <- mean(data_refe_10[!is.na(return_control_data$data[, i_series]) & !is.na(data_refe_10)])

              data_refe_ok <- (mean_10 / mean_refe_10) * data_refe_10
              return_control_data$data[is.na(return_control_data$data[, i_series]), i_series] <- data_refe_ok[is.na(return_control_data$data[, i_series])]
            }
          }
        }
      }
    }
  }
  return(return_control_data)
}

#' Rellena la serie recibida utilizando las otras en el orden en el que están en other_series
#'
#' @param series serie de datos a completar 
#' @param other_series series de datos con las que completar en el orden en el que se tienen que utilizar 
#'
#' @return serie de datosd rellena
#' @export
#'
fill_one_series <- function(series, other_series){

  return_series <- series

  # Por cercaniá, solo miramos correlación para eliminar, 5 años comunes

  i_other_series = 1
  while(sum(is.na(return_series)) > 0 & i_other_series <= dim(other_series)[2]){

    use_series <- other_series[, i_other_series]

    # Datos NA que son 0 en la serie de relleno se ponen directamente a 0
    return_series[is.na(return_series) & !is.na(use_series) & use_series == 0] = 0

    use_series_no_0 <- use_series[(is.na(use_series) | use_series != 0) & (is.na(series) | series != 0)]
    series_no_0 <- series[(is.na(use_series) | use_series != 0) & (is.na(series) | series != 0)]

    use_series_comoon = use_series_no_0[!is.na(series_no_0) & !is.na(use_series_no_0)]
    series_comoon = series_no_0[!is.na(series_no_0) & !is.na(use_series_no_0)]

    # Calculamos sobre los datos en común de la serie a usar para rellenar
    l_mom <- suppressWarnings(lmom::samlmu(use_series_comoon, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0))

    # Seleccionamos la mejor de las 10 funciones
    p.value = array(NA, dim = length(pel_functions), dimnames = list(names(pel_functions)))
    best_funcion = names(pel_functions)[3]
    for(best_funcion in names(pel_functions)){
      shapiro_value <- tryCatch({
        par_exp <- suppressWarnings(pel_functions[[best_funcion]](l_mom))
        cdf_exp <- suppressWarnings(cdf_functions[[best_funcion]](use_series_comoon, para = par_exp))
        norm_exp <- suppressWarnings(stats::qnorm(cdf_exp))
        shapiro_value <- suppressWarnings(stats::shapiro.test(norm_exp)$p.value)
      }, error = function(cond) {
        NA
      })
      p.value[best_funcion] <- shapiro_value
    }
    if(sum(!is.na(p.value)) > 0){
      best_funcion <- names(pel_functions)[which(p.value == max(p.value, na.rm = TRUE))[1]]

      # Utilizamos la mejor de las 10 funciones
      par_exp <- pel_functions[[best_funcion]](l_mom)

      # Calculamos sobre todos los datos de la serie a usar para rellenar
      cdf_exp <- cdf_functions[[best_funcion]](use_series_no_0, para = par_exp)

      generate_series <- suppressWarnings(as.vector(stats::quantile(series_no_0, cdf_exp, na.rm = TRUE)))
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
#'
#' @return data y coor con los datos que no pasan el control eliminados
#' @export
#'
save_data <- function(data_ori, control_data){
  data_return <- list()

  i_years <- ceiling(0.75 * (i_end - i_inis + 1))
  names(i_years) <- i_inis

  i_ini <- i_inis[length(i_inis)]
  for(i_ini in i_inis){
    data_ori_select <- data_ori[, colnames(control_data$data), drop = FALSE]

    i_year <- i_years[as.character(i_ini)]
    date <- seq(chron::chron(paste0("01/01/", i_ini), format=c(dates = "d/m/y", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s")), chron::chron(paste0("31/12/", i_end), format=c(dates = "d/m/y", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s")), by = "month")

    data <- array(NA, dim = c(length(date), dim(control_data$data)[2]))
    colnames(data) <- colnames(control_data$data)
    rownames(data) <- as.character(date)
    date_ok <- c(rownames(data), rownames(data_ori_select))[duplicated(c(rownames(data), rownames(data_ori_select)))]
    data[date_ok, ] <- as.matrix(data_ori_select)[date_ok, ]
    data_save <- data[, apply(data, c(2), sum_no_nas) > i_year * 12, drop = FALSE]
    coor_save <- control_data$coor[colnames(data_save), , drop = FALSE]
    data_save[date_ok, colnames(data_save)] <- control_data$data[date_ok, colnames(data_save), drop = FALSE]
    if(dim(coor_save)[1] > 0){      
      data_return[[paste0("start_", i_ini)]] <- list(data_save = data_save, coor_save = coor_save)
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
    data_save_year_month["year"] = as.numeric(substr(rownames(data_save_year_month), 8, 11))
    data_save_year_month["month"] = as.numeric(months(chron::chron(rownames(data_save_year_month), format=c(dates = "d/m/yy", times = "h:m:s"), out.format=c(dates = "d/m/yy", times = "h:m:s"))))
    data_save_year_month = data_save_year_month[c("year", "month", colnames(data_save_year_month)[1:(dim(data_save_year_month)[2]-2)])]

    utils::write.table(round(data_save_year_month, digits = 3), file.path(folder_name, paste0("data_", i_ini, ".csv")), row.names = FALSE, sep = ";")

    coor_save_csv <- as.data.frame(coor_save)
    coor_save_csv[, "coor"] <- rownames(coor_save_csv)
    coor_save_csv <- coor_save_csv[, c("coor", "lat", "lon")]
    coor_save_csv[, c("lat", "lon")] <- round(coor_save_csv[, c("lat", "lon")], digits = 3)
    utils::write.table(coor_save_csv, file.path(folder_name, paste0("coor_", i_ini, ".csv")), sep = ";", row.names = FALSE)
  }
}

#' Leemos los datos desde los CSVs con el formato acordado
#' Los ficheros de entrada son 2 CSVs uno de coordenadas en grados (filas las estaciones y columnas lat y lon y otro de datos mensuales con fechas en filas y las estaciones en las columnas)
#'
#' @param file_data ruta del fichero de datos
#' @param file_coor ruta del fichero de coordenadas
#'
#' @return datos originales, datos de interes y coordenadas de las estaciones leidas
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

    data_ori[data_ori < 0] <- NA

    data_ori[, 1:length(data_ori)] <- sapply(data_ori[, 1:length(data_ori)], as.numeric)
    coor[, 1:length(coor)] <- sapply(coor[, 1:length(coor)], as.numeric)

    # Revisión de fechas
    dates <- seq(chron::chron(rownames(data_ori)[1], format = time_format, out.format = time_format), chron::chron(rownames(data_ori)[dim(data_ori)[1]], format = time_format, out.format = time_format), by = "month")
    data <- data_ori[as.character(dates), , drop = FALSE]

    # Revisión de que tenemos las coordenadas
    no_coor <- which(!colnames(data) %in% rownames(coor))
    if(length(no_coor) > 0){
      warning(paste("We eliminate the stations", paste(colnames(data)[no_coor], sep = " ", collapse = " "), "because the coordinates are not available.", sep = " ", collapse = " "))
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
  clevels <- c(90, 92, 94, 95, 97.5, 99)
  if(is.na(match(clevel, clevels))){
     stop("clevel must be 90, 92, 94, 95, 97.5 or 99")
  }
  n <- NROW(x)
  y1 <- x[, 2]
  Z <- scale(y1)[, 1]
  Tv <- c()

  for (v in 1:(n-1)) {
    z1 <- base::mean(Z[1:v], na.rm = T)
    z2 <- base::mean(Z[(v + 1):n], na.rm = T)
    Tvi <- (v * (z1^2)) + ((n - v) * (z2^2))
    Tv <- rbind(Tv, Tvi)
  }
  # Test Statistic (omit tail ends +-5 years)
  T0 <- max(Tv[5:(length(Tv) - 5)])
  T0x <- which(Tv == T0)
  Tc <- as.numeric(snht_cval[max(which(snht_cval[, 'n'] <= n)), paste('clevel', clevel, sep = "")])

  return(list(Tc = Tc, T0 = T0, T0x = T0x))
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
  mkTrend_pval <- mkTrend_data["corrected_p_value"]
  mkTrend_pvalue0 <- mkTrend_data["p_value"]
  mkTrend_pval[is.na(mkTrend_pval) | is.null(mkTrend_pval)] <- mkTrend_pvalue0[is.na(mkTrend_pval) | is.null(mkTrend_pval)]

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
  year <- read_years(rownames(data))
  data_year <- as.data.frame(data) 
  data_year <- stats::aggregate(data_year, by = list("year" = year), FUN = sum)
  rownames(data_year) <- paste0("year_", as.character(unique(year)))
  data_year$year <- NULL
  return(as.matrix(data_year))
}

#' Devuelve el slope z por años y estaciones
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
  sens_slope[["summer"]] <- apply(calc_data_year(data[months %in% c("Jun", "Jul", "Aug"), , drop = FALSE]), c(2), calc_function)
  sens_slope[["autumn"]] <- apply(calc_data_year(data[months %in% c("Sep", "Oct", "Nov"), , drop = FALSE]), c(2), calc_function)

  data_winter <- data[months %in% c("Dec", "Jan", "Feb"), , drop = FALSE]
  if(!grepl("Dec", rownames(data_winter)[1])){
    data_winter <- data_winter[3:dim(data_winter)[1], , drop = FALSE]
  }
  i_dec <- which(grepl("Dec", rownames(data_winter)))
  rownames(data_winter)[i_dec[1:(length(i_dec) - 1)]] <- rownames(data_winter)[i_dec[2:length(i_dec)]]
  if(grepl("Dec", rownames(data_winter)[dim(data_winter)[1]])){
    data_winter <- data_winter[1:(dim(data_winter)[1] - 1), , drop = FALSE]
  }

  sens_slope[["winter"]]  <- apply(calc_data_year(data_winter), c(2), calc_function)
  sens_slope[["spring"]] <- apply(calc_data_year(data[months %in% c("Mar", "Apr", "May"), , drop = FALSE]), c(2), calc_function)
  return(sens_slope)
}

#' Calcula pval (a veces no da resultado por temas de iteración) entonces coger pval0.
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
   if(sum(!is.na(ro[i])) > 0 && (ro[i] > sig || ro[i] < -sig)) {
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
   z <- suppressWarnings((S - 1) / sqrt(VS))
   z0 <- suppressWarnings((S - 1)/sqrt(var.S))
  } else {
   z <- suppressWarnings((S + 1) / sqrt(VS))
   z0 <- suppressWarnings((S + 1) / sqrt(var.S))
  }

  if(is.finite(z)){
    pval <- 2 * stats::pnorm(- abs(z))
    pval0 <- 2 * stats::pnorm(- abs(z0))
  }else{
    pval <- pval0 <- 1
  }
  p_values <- c(pval0, pval)
  names(p_values) <- c("p_value", "corrected_p_value")

  return(p_values)
}

#' regresión lineal de los datos contra los años
#'
#' @param data index
#'
#' @return lm
#' @export
#'
calc_mkTrend_slp <- function(data){
  years <- read_years(names(data))
  if(sum(data) != 0){
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
    years <- read_years(names(index))

    n <- length(index)
    below <- which(index < threshold)
    above <- which(index >= threshold)
    dry <- index
    dry[below] <- 1
    dry[above] <- NA
    index2 <- - dry * index
    duration1 <- rep(NA, n)
    magnitud1 <- rep(NA, n)
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
    if(length(cases) > 0){
      spells <- cbind(years[cases], duration2[cases], magnitud2[cases])
      aggr1 <- stats::aggregate(spells[, 2], list(spells[, 1]), sum)
      aggr2 <- stats::aggregate(spells[, 3], list(spells[, 1]), sum)
      dry_periods <- as.data.frame(cbind(aggr1, aggr2[, 2]))
      years_complete <- data.frame(unique(years))
      colnames(years_complete) <- c("Group.1")
      base_datos <- merge(years_complete, dry_periods, by = 'Group.1', all = TRUE)
      base_datos[is.na(base_datos)] <- 0
      base_datos[, 2] <- base_datos[, 2] + 1
      base_datos[, 3] <- base_datos[, 3] + 1
      series_dry <- base_datos
      colnames(series_dry) <- c("years","duration","magnitude")
      a.1 <- as.vector(Kendall::Kendall(series_dry[, 1], series_dry[, 2])$sl)
      a.2 <- as.vector(Kendall::Kendall(series_dry[, 1], series_dry[, 3])$sl)
      lm.1 <- stats::lm(series_dry[, 2] ~ series_dry[,1])
      kk.1 <- as.data.frame(stats::coef(lm.1))
      a.3 <- kk.1[2, ]
      lm.2 <- stats::lm(series_dry[, 3] ~ series_dry[,1])
      kk.2 <- as.data.frame(stats::coef(lm.2))
      a.4 <- kk.2[2, ] * 10
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
  years <- read_years(names(datos))

  if(length(years) != length(unique(years))){
    stop("Repeating years in mobile_trends")
  }

  matriz_p <- NA
  matriz_s <- NA
  matriz_per <- NA
  if(length(years) > 30){
    matriz_p <- matrix(NA, length(years) - 30 + 1, length(years) - 30 + 1)
    matriz_s <- matrix(NA, length(years) - 30 + 1, length(years) - 30 + 1)
    matriz_per <- matrix(NA, length(years) - 30 + 1, length(years) - 30 + 1)

    bb <- 1
    for (bb in c(1:dim(matriz_p)[1])){
      aa <- length(years) - bb + 1
      rr <- aa - 30 + 1
      for (i in 1:rr){
        if(sum(datos[i:aa]) != 0){
          matriz_p[bb, i] <- calc_mkTrend_pval(datos[i:aa])
          theil_sen_regression_datos <- datos[i:aa]
          theil_sen_regression_years <- years[i:aa]
          matriz_s[bb, i] <- RobustLinearReg::theil_sen_regression(theil_sen_regression_datos ~ theil_sen_regression_years)$coefficients[2]
          matriz_per[bb, i] <- calc_percentage(datos[i:aa], years = years[i:aa])
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


#' Coeficientes de variación, desviación estándar
#' https://fhernanb.github.io/Manual-de-R/varia.html
#'
#' @param x datos
#' @param na.rm Ignorara NAs
#'
#' @return percentage
#' @export
#'
coef_var <- function(x, na.rm = FALSE) {
  stats::sd(x, na.rm = na.rm) / base::mean(x, na.rm = na.rm)
}

#' Diferencia en porcentaje
#'
#' @param datos datos
#' @param years años
#'
#' @return percentage
#' @export
#'
calc_percentage <- function(datos, years = NA) {
  if(sum(!is.na(years)) <= 0){
    years <- read_years(names(datos))
  }
  intercept_t <- RobustLinearReg::theil_sen_regression(datos ~ years)$coefficients
  final <- years[length(years)] * intercept_t[2] + intercept_t[1]
  inicial <- years[1] * intercept_t[2] + intercept_t[1]
  percentage <- 100 * final / inicial - 100
  if(is.na(percentage)){
    percentage <- 0
  }
  return(percentage)
}

#' Devuelve el procentage de datos válidos que son 0s
#'
#' @param data datos
#'
#' @return percentage
#' @export
#'
percentage_of_zeros <- function(data) {
  return(100 * sum(data == 0, na.rm = TRUE) / sum(!is.na(data), na.rm = TRUE))
}

#' Elimina datos si tenemos 8 meses o más seguidos de 0s, si uno de los meses implicados tiene menos del 70 por ciento de ceros se elimina su dato. Elimina también datos si tenemos 5 meses o más seguidos de 0s, si todos los meses implicados tiene menos del 70 por ciento de ceros se eliminan todos.
#'
#' @param data datos
#'
#' @return datos con los grupos de 0s eliminados
#' @export
#'
delete_zero <- function(data) {

  dates <- chron::chron(rownames(data), format = time_format, out.format = time_format)

  # Porcentaje de 0s para cada mes y estación
  months_percentage <- array(NA, dim = c(dim(data)[2], 12))
  rownames(months_percentage) = colnames(data)

  month <- 1
  for(month in c(1:12)){
    ok_month <- as.numeric(base::months(dates)) == month
    data_month <- data[ok_month, , drop = FALSE]
    months_percentage[, month] <- apply(data_month, c(2), percentage_of_zeros)
  }

  data_zero <- data
  
  stations <- colnames(data_zero)
  station <- stations[1]
  for(station in stations){
    station_rle <- rle(data_zero[, station])
    position <- cumsum(station_rle$lengths)

    # Elimina datos si tenemos 8 meses o más seguidos de 0s, si uno de los meses implicados tiene menos del 70 por ciento de ceros se elimina su dato.
    zero_groups <- which(!is.na(station_rle$values) & station_rle$values == 0 & station_rle$lengths >= 8)
    zero_group <- zero_groups[1]
    for(zero_group in zero_groups){
        range <- c((position[zero_group] - station_rle$lengths[zero_group] + 1) : position[zero_group])
        data_zero[range[months_percentage[station, as.numeric(base::months(dates[range]))] < 70], station] <- NA
    }

    # Elimina datos si tenemos 5 meses o más seguidos de 0s, si todos los meses implicados tiene menos del 70 por ciento de ceros se eliminan todos.
    zero_groups <- which(!is.na(station_rle$values) & station_rle$values == 0 & station_rle$lengths >= 5 & station_rle$lengths < 8)
    zero_group <- zero_groups[1]
    for(zero_group in zero_groups){
        range <- c((position[zero_group] - station_rle$lengths[zero_group] + 1) : position[zero_group])
        if(length(range) == sum(months_percentage[station, as.numeric(base::months(dates[range]))] < 70)){
          data_zero[range[months_percentage[station, as.numeric(base::months(dates[range]))] < 70], station] <- NA
        }
    }

  }

  return(data_zero)
}

#'Para cada estación guardar, número de datos de entrada y de datos eliminados
#'
#' @param ori_data datos iniciales
#' @param process_data datos procesados
#' @param folder carpeta en la que se guarda el fichero resultante
#'
#' @return datos eliminados y datos de entrada por estación
#' @export
#'
save_delete_data <- function(ori_data, process_data, folder) {

  ori_data <- ori_data[rownames(process_data), colnames(process_data), drop = FALSE]

  ori_data[!is.na(ori_data)] <- 1
  ori_data[is.na(ori_data)] <- 0
  process_data[!is.na(process_data)] <- 1
  process_data[is.na(process_data)] <- 0
  sum_ori <- apply(ori_data, c(2), sum)
  sum_process <- apply(process_data, c(2), sum)
  change <- sum_ori - sum_process
  data_save <- array(NA, dim = c(length(change), 2))
  rownames(data_save) <- colnames(ori_data)
  colnames(data_save) <- c("data_ori", "data_delete")
  data_save[, "data_ori"] <- as.numeric(sum_ori)
  data_save[, "data_delete"] <- as.numeric(change)

  # dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  # utils::write.table(data_save, file.path(folder, "delete_data.csv"), sep = ";")

  return(data_save)
}

#'Calcular estadísticos de la reconstrucción
#' - hydroGOF -- estadísitico por estación
#' - D / MAE / PBIAS / RMSE  - por estación y mes
#'
#' @param sim datos rellenados
#' @param obs datos iniciales
#'
#' @return datos eliminados y datos de entrada por estación
#' @export
#'
calculate_reconstruction_statistics <- function(sim, obs) {
  reconstruction_statistics <- c()  
  reconstruction_statistics["d"] <- hydroGOF::d(sim, obs)
  reconstruction_statistics["mae"] <- hydroGOF::mae(sim, obs)
  reconstruction_statistics["pbias"] <- hydroGOF::pbias(sim, obs)
  reconstruction_statistics["rmse"] <- hydroGOF::rmse(sim, obs)
  return(reconstruction_statistics)
}
