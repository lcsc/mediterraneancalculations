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

# library(rworldmap)
# library(rworldxtra)

# library(testthat)
# library(mediterraneancalculations)

context("MediterraneanCalculations")

##################################################

max_dist <- 200 # Máxima distancia entre 2 estaciones para usar una para evaluar o completar la otra

#' genera un data asociado a un coor que sirve en varios tests
#'
#' @param i_year año de comienzo de los datos
#' @param ndata número de estaciones
#'
#' @return list de data y coor
#' @export
#'
generate_mock_data_coor <- function(i_year = 2001, ndata = 5){
  times <- as.character(seq(chron::chron(paste0("01/Jan/", i_year), format = time_format, out.format = time_format), chron::chron("21/Dec/2021", format = time_format, out.format = time_format), by = "month"))
  data <- array(1000, dim = c(length(times), ndata))
  colnames(data) <- paste0("X", c(1:dim(data)[2]))
  rownames(data) <- times
  data[12, 1] <- NA
  data[18, 1] <- NA
  data[19, 1] <- 3 # Dato a eliminar por el control
  data[20, 1] <- 2 
  # data[, 2] <- NA
  data[1:(8*12), 2] <- rep(2, 8*12) # Serie a eliminar por pocos datos
  data[, 3] <- rep(2, dim(data)[1])
  data[20, 3] <- 130 # Dato a eliminar por el control
  data[20, 4] <- 130
  data[20, 5] <- 130

  coor <- array(NA, dim = c(dim(data)[2], 2))
  rownames(coor) <- colnames(data)
  colnames(coor) <- c("lat", "lon")

  coor[1, ] <- c(1, 1)
  coor[2, ] <- c(1, 1.1)
  coor[3, ] <- c(1.1, 1)
  coor[4, ] <- c(1.1, 1.1)
  coor[5, ] <- c(1, 0.9)

  # Estación solitaria
  if(ndata > 5){
    data[, 6] <- data[, 5]
    coor[6, ] <- c(1.2, 1.2)
  }  
  if(ndata > 6){  
    data[, 7] <- data[, 5]
    coor[7, ] <- c(1.3, 1.3)
  }
  if(ndata > 7){
    data[, 8] <- data[, 5]
    coor[8, ] <- c(90, 45)
  }

  return(list(data = data, coor = coor))
}

# ############################################################

#' testea la función order_data
#'
#' @return None
#' @export
#'
test_that("order_data", {
  control_data <- generate_mock_data_coor()
  data_result <- order_data(t(control_data$data))
  expect_equivalent(c(20, 1), data_result$X5[c(1:2)], info = deparse(sys.calls()[[sys.nframe()]]))
})

#' testea la función sum_no_nas
#'
#' @return None
#' @export
#'
test_that("sum_no_nas", {
  data <- c(1:300)
  data[200:299] <- NA
  data_result <- sum_no_nas(data)
  expect_equivalent(200, data_result, info = deparse(sys.calls()[[sys.nframe()]]))
})

#' testea la función select_data
#'
#' @return None
#' @export
#'
test_that("select_data", {
  data <- c(1, NA, 2, NA, NA, 3, 4, 5, 6, 7, NA, 8)
  data_result <- select_data(data, 5)
  expect_equivalent(c(1:5), data_result, info = deparse(sys.calls()[[sys.nframe()]]))
})

#' testea la función apply_ecdf_month
#'
#' @return None
#' @export
#'
test_that("apply_ecdf_month", {
  control_data <- generate_mock_data_coor()
  control_data$data[13, ] = c(34, 34, 34, 34, 34)
  data_result <- apply_ecdf_month(control_data$data[grepl("Jan", rownames(control_data$data)), ])
  expect_equivalent(data_result[2, "X4"] < data_result[3, "X4"], TRUE)
})

#' testea la función apply_ecdf
#'
#' @return None
#' @export
#'
test_that("apply_ecdf", {
  control_data <- generate_mock_data_coor()
  control_data$data[13, ] = c(34, 34, 34, 34, 34)
  data_result <- apply_ecdf(control_data$data)
  expect_equivalent(data_result[13, "X4"] < data_result[14, "X4"], TRUE)
})

#' testea la función near_correlations
#'
#' @return None
#' @export
#'
test_that("near_correlations", {
  control_data <- generate_mock_data_coor()

  control_data$data[11, 1] <- 3
  control_data$data[11, 2] <- -3
  control_data$data[11, 3] <- -1
  control_data$data[11, 4] <- 1
  control_data$data[11, 5] <- 0
  # print(near_correlations)
  data_result <- near_correlations(control_data$data, control_data$coor, max_dist = max_dist)
  expect_equivalent(is.na(NA), is.na(data_result$Aug[3, 3]), info = deparse(sys.calls()[[sys.nframe()]]))
  expect_equivalent(-1, data_result$Aug[1, 3], info = deparse(sys.calls()[[sys.nframe()]]))
})

#' testea la función near_estations
#'
#' @return None
#' @export
#'
test_that("near_estations", {
  control_data <- generate_mock_data_coor()
  control_data$coor[5, ] <- c(50, 50)
  data_result <- near_estations(data = control_data$data, coor = control_data$coor, max_dist = max_dist)
  expect_equivalent(0, length(data_result$X5), info = deparse(sys.calls()[[sys.nframe()]]))
  expect_equivalent(3, length(data_result$X1), info = deparse(sys.calls()[[sys.nframe()]]))

  control_data$coor[1, ] <- c(1, 1)
  control_data$coor[2, ] <- c(2, 1)
  control_data$coor[3, ] <- c(3, 1)
  control_data$coor[4, ] <- c(4, 1)
  control_data$coor[5, ] <- c(5, 1)
  data_result <- near_estations(data = control_data$data, coor = control_data$coor, max_dist = max_dist)
  expect_equivalent(4, length(data_result$X5), info = deparse(sys.calls()[[sys.nframe()]]))
  expect_equivalent(c(2, 3, 4, 5), data_result$X1, info = deparse(sys.calls()[[sys.nframe()]]))
})

#' testea la función quality_control
#'
#' @return None
#' @export
#'
test_that("quality_control", {
  control_data <- generate_mock_data_coor()
  control_data$data[, 2] <- NA

  control_data$data[31, 1] <- 3
  control_data$data[31, 2] <- -3
  control_data$data[31, 3] <- -1
  control_data$data[31, 4] <- 1
  control_data$data[31, 5] <- 0

  control_data$data[32, 5] <- 0

  data_result <- control_data$data[, c(1, 3, 4, 5)]
  data_result[19, 1] <- NA
  data_result[20, 2] <- NA
  data_result[32, 4] <- NA

  coor_result <- control_data$coor[colnames(data_result), ]

  result <- quality_control(data = control_data$data, coor = control_data$coor, max_dist = max_dist)

  expect_equivalent(result$data, data_result, info = deparse(sys.calls()[[sys.nframe()]]))
  expect_equivalent(result$coor, coor_result, info = deparse(sys.calls()[[sys.nframe()]]))

  # Estación solitaria
  control_data$data[12, "X5"] <- NA
  control_data$coor["X5", ] <- c(45, 90)
  result <- quality_control(data = control_data$data, coor = control_data$coor, max_dist = max_dist)

  expect_equivalent(result$data[, "X5"], control_data$data[, "X5"], info = deparse(sys.calls()[[sys.nframe()]]))
})

#' testea la función overlap_station
#'
#' @return None
#' @export
#'
test_that("overlap_station", {
  control_data <- generate_mock_data_coor()

  control_data$data[, ] <- 1000
  control_data$data[18, 2] <- NA
  control_data$data[19, 2] <- NA

  data_ok <- array(NA, dim = c(dim(control_data$data)[2], dim(control_data$data)[2]), dimnames = list(colnames(control_data$data), colnames(control_data$data)))
  data_ok[, ] <- dim(control_data$data)[1]
  data_ok[, 2] <- dim(control_data$data)[1] - 2 
  data_ok[2, ] <- dim(control_data$data)[1] - 2 
  data_result <- overlap_station(control_data = list(data = control_data$data, coor = control_data$coor))
  expect_equivalent(data_result, data_ok, info = deparse(sys.calls()[[sys.nframe()]]))
})

#' testea la función fill_one_series
#'
#' @return None
#' @export
#'
test_that("fill_one_series", {
  info = deparse(sys.calls()[[sys.nframe()]])
  control_data <- generate_mock_data_coor()

  control_data$data[12, 1] <- NA
  control_data$data[24, 2] <- NA
  control_data$data[1, 3] <- -NA
  control_data$data[30, 4] <- NA
  control_data$data[60, 5] <- NA

  data_result <- fill_one_series(series = control_data$data[, 1], other_series = control_data$data[, -1])
  expect_equivalent(sum(is.na(data_result)), 0, info = info)
})

#' testea la función fill_series
#'
#' @return None
#' @export
#'
test_that("fill_series", {
  info = deparse(sys.calls()[[sys.nframe()]])
  control_data <- generate_mock_data_coor()

  control_data$data[24, 1] <- 3
  control_data$data[24, 2] <- -3
  control_data$data[24, 3] <- -1
  control_data$data[24, 4] <- 1
  control_data$data[24, 5] <- 0

  data_result <- fill_series(control_data, min_correlation = min_correlation, max_dist = max_dist)
  expect_lt(sum(is.na(data_result$data)), sum(is.na(control_data$data)), label = info)
  expect_equivalent(sum(is.na(data_result$data)), 0, info = info)
})

#' testea la función save_data
#'
#' @return None
#' @export
#'
test_that("save_data", {
  control_data <- generate_mock_data_coor(i_year = 1950) 
  data_ori <- control_data$data
  data_ori[1, 1] <- NA
  data_ori[12, 4] <- NA

  data <- save_data(data_ori, control_data)

  expect_equivalent(dim(data[["start_1930"]]$data)[2], 5)
  expect_equivalent(dim(data[["start_1930"]]$coor)[1], 5)
})

#' testea la función save_csvs and read_data
#'
#' @return None
#' @export
#'
test_that("save_csvs_read_data", {
  control_data <- generate_mock_data_coor(i_year = 1950) 

  i_ini <- "1930"
  folder_name <- "test_result"
  save_csvs(i_ini, folder_name, data_save = control_data$data, coor_save = control_data$coor)
  data <- read_data(file_data = "test_result/data_1930.csv", file_coor = "test_result/coor_1930.csv")

  expect_equivalent(as.data.frame(control_data$data), data$data)
  expect_equivalent(as.data.frame(control_data$coor), data$coor)
})

#' testea la función alexanderson_homogenize_data
#'
#' @return None
#' @export
#'
test_that("alexanderson_homogenize_data", {
  file_data <- generate_mock_data_coor()
  file_data$data[c(1:100), 1] <- c(1:100)
  file_data$data[c(101:dim(file_data$data)[1]), 1] <- c(-101:-dim(file_data$data)[1])
  data <- alexanderson_homogenize_data(file_data)
  expect_lt(data[1, 1], 0)

  # Estación solitaria
  file_data$data[12, "X5"] <- NA
  file_data$coor["X5", ] <- c(45, 90)
  data <- alexanderson_homogenize_data(file_data)
  expect_lt(data[1, 1], 0)
})

#' testea la función save_parameters
#'
#' @return None
#' @export
#'
test_that("save_parameters", {
  file_data <- generate_mock_data_coor()
  if(file.exists(file.path("test_result", "parameters.txt"))){
    file.remove(file.path("test_result", "parameters.txt"))
  }
  save_parameters("test_result", 200)
  expect_equivalent(file.exists(file.path("test_result", "parameters.txt")), TRUE)
})

#' testea la función second_data_fill_data
#'
#' @return None
#' @export
#'
test_that("second_data_fill_data", {
  file_data <- generate_mock_data_coor()
  file_data$data[14, "X2"] <- NA
  file_data$data[140, "X3"] <- NA
  file_data$data[100, "X4"] <- NA
  file_data$data[200, "X5"] <- NA # Aug
  data <- second_data_fill_data(file_data)
  expect_equivalent(file_data$coor[c("X3", "X5"), ], data$coor)

  # Estación solitaria
  # file_data$data[12, "X5"] <- NA
  file_data$coor["X5", ] <- c(45, 90)
  data <- second_data_fill_data(file_data)
  expect_equivalent(file_data$coor[c("X3", "X5"), ], data$coor)
})

#' testea la función calculate_statistics_data
#'
#' @return None
#' @export
#'
test_that("calculate_statistics_data", {
  file_data <- generate_mock_data_coor()
  file_data$data[is.na(file_data$data)] <- 35
  data <- calculate_statistics_data(file_data)
  expect_equivalent(data$average[12], mean(file_data$data[12, ], na.rm = TRUE))
})

#' testea la función dry_spell_trend
#'
#' @return None
#' @export
#'
test_that("dry_spell_trend", {
  file_data <- generate_mock_data_coor(i_year = 1950)
  file_data$data[is.na(file_data$data)] <- 35
  set.seed(42)
  file_data$data[, 3] <- sample(c(dim(file_data$data)[1]:1))/400 - 2
  data <- dry_spell_trend(index = file_data$data[, 3], threshold = 0)
  # expect_equivalent(round(data, 5), c(0.67246, 0.88021, 0, -7e-05)) # Función de Sergio, comprobamos que mantiene resultado
  expect_equivalent(round(data, 3), c(0.265, 0.137, 0, 0)) # Función de Sergio, comprobamos que mantiene resultado
})

#' testea la función mobile_trends
#'
#' @return None
#' @export
#'
test_that("mobile_trends", {
  file_data <- generate_mock_data_coor(i_year = 1950)
  file_data$data[is.na(file_data$data)] <- 35
  data <- mobile_trends(datos = file_data$data[, 3])
  expect_equivalent(dim(data$matriz_s)[1], 43)
  expect_equivalent(round(mean(data$matriz_p, na.rm = TRUE), 5), 0.29986) # Función de Sergio, comprobamos que mantiene resultado
})

#' testea la función mkTrend
#'
#' @return None
#' @export
#'
test_that("mkTrend", {
  file_data <- generate_mock_data_coor(i_year = 1950)
  data <- mkTrend(x = file_data$data[, 3])
  expect_equivalent(round(mean(as.numeric(as.matrix(data)), na.rm = TRUE), 5), 27.8223) # Función de Sergio, comprobamos que mantiene resultado
})

#' testea la función calc_data_year_month_station
#'
#' @return None
#' @export
#'
test_that("calc_data_year_month_station", {
  file_data <- generate_mock_data_coor()
  data <- unlist(calc_data_year_month_station(data = file_data$data, mean))
  expect_equivalent(names(which(is.na(data))), c("year.X1", "summer.X1", "winter.X1", "Jun.X1", "Dec.X1"))
})

#' testea la función calc_mkTrend_pval
#'
#' @return None
#' @export
#'
test_that("calc_mkTrend_pval", {
  file_data <- generate_mock_data_coor()
  data <- calc_mkTrend_pval(data = file_data$data[, 2])
  expect_equivalent(round(data, 5), 0.00131) # Función de Sergio, comprobamos que mantiene resultado
})

#' testea la función snht
#'
#' @return None
#' @export
#'
test_that("snht", {
  file_data <- generate_mock_data_coor()
  data <- snht(x = file_data$data, clevel = 95)$T0x
  expect_equivalent(data, 96) # Función de Sergio, comprobamos que mantiene resultado
})

#' testea la función calc_data_year
#'
#' @return None
#' @export
#'
test_that("calc_data_year", {
  file_data <- generate_mock_data_coor()
  data <- calc_data_year(file_data$data)
  expect_equivalent(data["year_2001", "X3"], sum(file_data$data[grepl("2001", rownames(file_data$data)), "X3"])) 
})

#' testea la función fill_unfillable_station
#'
#' @return None
#' @export
#'
test_that("fill_unfillable_station", {
  file_data <- generate_mock_data_coor(i_year = 1980)
  fillable_years <- 4
  file_data$data[100, 1] <- NA
  file_data$data[140, 1] <- NA
  file_data$data[160, 1] <- NA
  file_data$data[1, 2] <- NA
  file_data$data[140, 3] <- NA
  file_data$data[160, 3] <- NA
  file_data$data[140, 4] <- NA
  file_data$data[160, 5] <- NA
  data <- fill_unfillable_station(data = file_data, fillable_years = fillable_years)
  expect_equivalent(sum(is.na(file_data$data[, 1])), sum(is.na(data$data[, 1])))
  expect_equivalent(sum(is.na(file_data$data[, 2])), sum(is.na(data$data[, 2])))
  expect_equivalent(sum(is.na(data$data[, 3])), 0) 
  expect_equivalent(sum(is.na(data$data[, 4])), 0) 
  expect_equivalent(sum(is.na(data$data[, 5])), 0) 
})

# #' testea la función main_mediterranean_calculations
# #'
# #' @return None
# #' @export
# #'
# test_that("main_mediterranean_calculations", {
#   control_data <- generate_mock_data_coor(i_year = 1950, ndata = 8) 
#   i_ini <- "1930"
#   folder_name <- "test_result"
#   save_csvs(i_ini, folder_name, data_save = control_data$data, coor_save = control_data$coor)

#   file_data <- "test_result/data_1930.csv"
#   file_coor <- "test_result/coor_1930.csv"

#   main_mediterranean_calculations(file_data = file_data, file_coor = file_coor, max_dist = max_dist)

#   # sum(is.na(control_data$data[,"X5"]))

# })

# #' testea la función main_mediterranean_calculations con 6 datasets
# #'
# #' @return None
# #' @export
# #'
# test_that("main_mediterranean_calculations", {
#   control_data <- generate_mock_data_coor(i_year = 1950, ndata = 6) 
#   i_ini <- "1930"
#   folder_name <- "test_result"
#   save_csvs(i_ini, folder_name, data_save = control_data$data, coor_save = control_data$coor)

#   file_data <- "test_result/data_1930.csv"
#   file_coor <- "test_result/coor_1930.csv"

#   main_mediterranean_calculations(file_data = file_data, file_coor = file_coor, max_dist = max_dist)

#   # sum(is.na(control_data$data[,"X5"]))

# })

# calc_mkTrend_slp

# #' tes pruebas
# #'
# #' @return None
# #' @export
# #'
# prueba_data <- function(){
#   file_data <- file.path("israel", "israel_summer_data.csv")
#   data_ori <- read.table(file_data, sep = ";", header = TRUE)
#   data <- read.table("/mnt/dostb1/DATOS/sergio/calculos_mediterraneo/results/data_1980.csv", sep = ";", header = TRUE)

#   rownames(data_ori) <- as.character(chron(paste0("01", "/", data_ori[, "month"], "/", data_ori[, "year"]), format = c(dates = "d/m/y", times = "h:m:s"), out.format = time_format))
#   data_ori["month"] <- NULL
#   data_ori["year"] <- NULL

#   rownames(data) <- as.character(chron(paste0("01", "/", data[, "month"], "/", data[, "year"]), format = c(dates = "d/m/y", times = "h:m:s"), out.format = time_format))
#   data["month"] <- NULL
#   data["year"] <- NULL


#   data_ori <- data_ori[rownames(data), colnames(data)]

#   station <- colnames(data_ori)[129]
#   sum(data[, station] != data_ori[, station], na.rm = TRUE)
#   sum(is.na(data[, station]))
#   sum(is.na(data_ori[, station])) 

#   # which(!is.na(data_ori["01/04/1986", ]) & !is.na(data_ori["01/10/1984", ]) & is.na(data["01/04/1986", ]) & is.na(data["01/10/1984", ]))
# }

