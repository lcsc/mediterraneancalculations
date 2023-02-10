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

max_dist_test <- 200 # Máxima distancia entre 2 estaciones para usar una para evaluar o completar la otra

#' genera un data asociado a un coor que sirve en varios tests
#'
#' @param i_year año de comienzo de los datos
#' @param ndata número de estaciones
#'
#' @return list de data y coor
#' @export
#'
generate_mock_data_coor <- function(i_year = 2001, ndata = 5){
  times <- as.character(seq(chron::chron(paste0("01/Jan/", i_year), format = time_format, out.format = time_format), chron::chron(paste0("21/Dec/", i_end), format = time_format, out.format = time_format), by = "month"))
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

  coor <- array(NA, dim = c(dim(data)[2], 2))
  rownames(coor) <- colnames(data)
  colnames(coor) <- c("lat", "lon")

  coor[1, ] <- c(1, 1)
  coor[2, ] <- c(1, 1.1)
  coor[3, ] <- c(1.1, 1)

  if(ndata > 3){    
    data[20, 3] <- 130 # Dato a eliminar por el control
    data[20, 4] <- 130
    data[20, 5] <- 130
    
    coor[4, ] <- c(1.1, 1.1)
    coor[5, ] <- c(1, 0.9)
  }

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

#' testea la función read_years
#'
#' @return None
#' @export
#'
test_that("read_years", {
  file_data <- generate_mock_data_coor()
  data <- read_years(txt = rownames(file_data$data))
  expect_equivalent(data[123], 2011) 
})

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
  data_result <- near_correlations(control_data$data, control_data$coor, max_dist = max_dist_test)
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
  data_result <- near_estations(data = control_data$data, coor = control_data$coor, max_dist = max_dist_test)
  expect_equivalent(0, length(data_result$X5), info = deparse(sys.calls()[[sys.nframe()]]))
  expect_equivalent(3, length(data_result$X1), info = deparse(sys.calls()[[sys.nframe()]]))

  control_data$coor[1, ] <- c(1, 1)
  control_data$coor[2, ] <- c(2, 1)
  control_data$coor[3, ] <- c(3, 1)
  control_data$coor[4, ] <- c(4, 1)
  control_data$coor[5, ] <- c(5, 1)
  data_result <- near_estations(data = control_data$data, coor = control_data$coor, max_dist = max_dist_test)
  expect_equivalent(4, length(data_result$X5), info = deparse(sys.calls()[[sys.nframe()]]))
  expect_equivalent(c(2, 3, 4, 5), data_result$X1, info = deparse(sys.calls()[[sys.nframe()]]))
})

#' testea la función quality_control
#'
#' @return None
#' @export
#'
test_that("quality_control", {

  max_diff_anomaly <- 0.6 # Máxima diferencia de anomalías para mantener dato en el control
  max_diff_anomaly_0 <- 0.5 # Máxima diferencia de anomalías para mantener dato en el control, si el dato es 0

  control_data <- generate_mock_data_coor(ndata = 6)
  control_data$data[, 2] <- NA

  control_data$data[31, 1] <- 3
  control_data$data[31, 2] <- -3
  control_data$data[31, 3] <- -1
  control_data$data[31, 4] <- 1
  control_data$data[31, 5] <- 0

  control_data$data[32, 5] <- 0

  data_result <- control_data$data[, c(1, 3, 4, 5, 6)]
  # data_result[19, 1] <- NA # ¿?
  data_result[20, 2] <- NA
  data_result[32, 4] <- NA

  coor_result <- control_data$coor[colnames(data_result), ]

  result <- quality_control(data = control_data$data, coor = control_data$coor, max_dist = max_dist_test, max_diff_anomaly = max_diff_anomaly, max_diff_anomaly_0 = max_diff_anomaly_0)

  expect_equivalent(result$data, data_result, info = deparse(sys.calls()[[sys.nframe()]]))
  expect_equivalent(result$coor, coor_result, info = deparse(sys.calls()[[sys.nframe()]]))

  control_data <- generate_mock_data_coor(ndata = 6)
  control_data$data[, "X1"] <- 1000:(1000 + dim(control_data$data)[1] - 1)
  control_data$data[, "X2"] <- 1000:(1000 + dim(control_data$data)[1] - 1)
  control_data$data[, "X3"] <- 1000:(1000 + dim(control_data$data)[1] - 1)
  control_data$data[, "X4"] <- 1000:(1000 + dim(control_data$data)[1] - 1)
  control_data$data[, "X5"] <- 1000:(1000 + dim(control_data$data)[1] - 1)
  control_data$data[, "X6"] <- 1000:(1000 + dim(control_data$data)[1] - 1)
  control_data$data[dim(control_data$data)[1], "X2"] <- 1
  result <- quality_control(data = control_data$data, coor = control_data$coor, max_dist = max_dist_test, max_diff_anomaly = max_diff_anomaly, max_diff_anomaly_0 = max_diff_anomaly_0)

  # Estación solitaria
  control_data$data[12, "X5"] <- NA
  control_data$coor["X5", ] <- c(45, 90)
  result <- quality_control(data = control_data$data, coor = control_data$coor, max_dist = max_dist_test, max_diff_anomaly = max_diff_anomaly, max_diff_anomaly_0 = max_diff_anomaly_0)

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

  data_result <- fill_series(control_data = control_data, min_correlation = min_correlation, max_dist = max_dist_test)
  expect_lt(sum(is.na(data_result$data)), sum(is.na(control_data$data)), label = info)
  expect_equivalent(sum(is.na(data_result$data)), 0, info = info)
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

  expect_equivalent(dim(data[["start_1931"]]$data)[2], 5)
  expect_equivalent(dim(data[["start_1931"]]$coor)[1], 5)
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

#' testea la función calc_mkTrend_pval
#'
#' @return None
#' @export
#'
test_that("calc_mkTrend_pval", {
  file_data <- generate_mock_data_coor()
  ok_data <- calc_data_year(data = file_data$data)
  data <- calc_mkTrend_pval(data = ok_data[, 2])
  expect_equivalent(round(data, 5), 0.00614) # 0.01054 ¿? Función de Sergio, comprobamos que mantiene resultado

  file_data_one <- generate_mock_data_coor(i_year = 1981, ndata = 3) 
  ok_data <- calc_data_year(data = file_data$data)
  ok_data[, "X1"] <- 1:length(ok_data[, "X1"])
  ok_data[, "X2"] <- length(ok_data[, "X2"]):1
  ok_data[, "X3"] <- 1 

  data <- calc_mkTrend_pval(data = ok_data[, "X1"])
  expect_equivalent(round(data, 3), c(0)) 
  data <- calc_mkTrend_pval(data = ok_data[, "X2"])
  expect_equivalent(round(data, 3), c(0)) 
  data <- calc_mkTrend_pval(data = ok_data[, "X3"]) 
  expect_equivalent(round(data, 3), c(1))
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

#' testea la función calc_data_year_month_station
#'
#' @return None
#' @export
#'
test_that("calc_data_year_month_station", {
  file_data <- generate_mock_data_coor()
  data <- unlist(calc_data_year_month_station(data = file_data$data, calc_function = base::mean))
  expect_equivalent(names(which(is.na(data))), c("year.X1", "summer.X1", "winter.X1"))

  file_data <- generate_mock_data_coor()
  file_data$data[c("01/Jun/2015", "01/Jul/2015", "01/Aug/2015"), ] <- 2
  file_data$data[c("01/Dec/2019", "01/Jan/2020", "01/Feb/2020"), ] <- 1
  data <- unlist(calc_data_year_month_station(data = file_data$data, calc_function = base::min))
  expect_equivalent(names(which(is.na(data))), c("year.X1", "summer.X1", "winter.X1"))
})

#' testea la función mkTrend
#'
#' @return None
#' @export
#'
test_that("mkTrend", {
  # Función original mkTrend
  mkTrend_ori <- function(x, ci = .95) {
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

  file_data <- generate_mock_data_coor(i_year = 1950)
  data <- mkTrend(x = file_data$data[, 3])
  data_pval <- calc_mkTrend_pval(data = file_data$data[, 3])
  data_ori <- mkTrend_ori(x = file_data$data[, 3])
  expect_equivalent(round(mean(as.numeric(as.matrix(data)), na.rm = TRUE), 3), 0.502)
  expect_equivalent(round(mean(as.numeric(as.matrix(data)), na.rm = TRUE), 5), round(mean(c(data_ori$p.value, data_ori$`Corrected p.value`), na.rm = TRUE), 5))  # Función de Sergio, comprobamos que mantiene resultado de su versión original
  expect_equivalent(round(data["corrected_p_value"], 5), round(data_pval, 5))
})

#' testea la función calc_mkTrend_slp
#'
#' @return None
#' @export
#'
test_that("calc_mkTrend_slp", {
  file_data <- generate_mock_data_coor()
  ok_data <- calc_data_year(data = file_data$data)
  data <- calc_mkTrend_slp(data = ok_data[, 2])
  expect_equivalent(round(data, 3), 647.825) # 0 ¿? Función de Sergio, comprobamos que mantiene resultado

  # Prueba pedida por Sergio
  file_data <- generate_mock_data_coor()

  file_data$data[, 1] <- c(1:dim(file_data$data)[1])
  file_data$data[, 2] <- rep(1, dim(file_data$data)[1])
  file_data$data[, 3] <- c(dim(file_data$data)[1]:1)
  ok_data <- calc_data_year(data = file_data$data)
  data1 <- calc_mkTrend_slp(data = ok_data[, 1])
  data2 <- calc_mkTrend_slp(data = ok_data[, 2])
  data3 <- calc_mkTrend_slp(data = ok_data[, 3])

  expect_false(data1 == 1)
  expect_true(data2 == 0)
  expect_false(data3 == -1)

  ok_data[, 1] <- c(1:dim(ok_data)[1])
  ok_data[, 2] <- rep(1, dim(ok_data)[1])
  ok_data[, 3] <- c(dim(ok_data)[1]:1)
  data1 <- calc_mkTrend_slp(data = ok_data[, 1])
  data2 <- calc_mkTrend_slp(data = ok_data[, 2])
  data3 <- calc_mkTrend_slp(data = ok_data[, 3])
  expect_true(data1 == 1)
  expect_true(data2 == 0)
  expect_true(data3 == -1)

  file_data_one <- generate_mock_data_coor(i_year = 1981, ndata = 3) 
  ok_data <- calc_data_year(data = file_data$data)
  ok_data[, "X1"] <- 1:length(ok_data[, "X1"])
  ok_data[, "X2"] <- length(ok_data[, "X2"]):1
  ok_data[, "X3"] <- 1

  data <- calc_mkTrend_slp(data = ok_data[, "X1"])
  expect_equivalent(round(data, 3), c(1)) 
  data <- calc_mkTrend_slp(data = ok_data[, "X2"])
  expect_equivalent(round(data, 3), c(-1)) 
  data <- calc_mkTrend_slp(data = ok_data[, "X3"])
  expect_equivalent(round(data, 3), c(0)) 
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
  # expect_equivalent(round(data, 3), c(0.265, 0.137, 0, 0)) # Función de Sergio, comprobamos que mantiene resultado
  expect_equivalent(round(data, 2), c(0.34, 0.28, 0.03, 0.45)) # Función de Sergio, comprobamos que mantiene resultado
})

#' testea la función mobile_trends
#'
#' @return None
#' @export
#'
test_that("mobile_trends", {
  file_data <- generate_mock_data_coor(i_year = 1950)
  file_data$data[is.na(file_data$data)] <- 35
  # system.time(mobile_trends(datos = file_data$data[, 3])) # 8.879   0.004   8.893
  data <- mobile_trends(datos = file_data$data[seq(1, dim(file_data$data)[1], by = 12), 3])
  expect_equivalent(dim(data$matriz_s)[1], 42) # Antes era 43
  expect_equivalent(round(mean(data$matriz_p, na.rm = TRUE), 5), 1) # Función de Sergio, comprobamos que mantiene resultado

  file_data <- generate_mock_data_coor(i_year = 1981, ndata = 3) # Secalcula para c(1871, 1901, 1931) 
  ok_data <- calc_data_year(data = file_data$data)
  ok_data[, "X1"] <- 1:length(ok_data[, "X1"])
  ok_data[, "X2"] <- length(ok_data[, "X2"]):1
  ok_data[, "X3"] <- 1
  data1 <- mobile_trends(datos = ok_data[, "X1"]) # 0s 0s 1s
  data2 <- mobile_trends(datos = ok_data[, "X2"]) # 0s 0s 1s
  data3 <- mobile_trends(datos = ok_data[, "X3"]) # 0s 0s 1s
  expect_equivalent(round(data1$matriz_p[1, 1], 3), 0)
  expect_equivalent(round(data2$matriz_p[1, 1], 3), 0)
  expect_equivalent(round(data3$matriz_p[1, 1], 3), 1)
})

#' testea la función calc_percentage
#'
#' @return None
#' @export
#'
test_that("calc_percentage", {
  file_data <- generate_mock_data_coor()
  ok_data <- calc_data_year(data = file_data$data)
  data <- calc_percentage(datos = ok_data[, "X4"])
  expect_equivalent(data, 0) # Función de Sergio, comprobamos que mantiene resultado

  file_data_one <- generate_mock_data_coor(i_year = 1981, ndata = 3) 
  ok_data <- calc_data_year(data = file_data_one$data)[, 1:3]
  ok_data[, "X1"] <- 1:length(ok_data[, "X1"])
  ok_data[, "X2"] <- length(ok_data[, "X2"]):1
  ok_data[, "X3"] <- 1

  data1 <- calc_percentage(datos = ok_data[, "X1"]) #MAL 35.13 -97.5 0 
  expect_equivalent(round(data1, 3), c(3900)) 
  data2 <- calc_percentage(datos = ok_data[, "X2"])
  expect_equivalent(round(data2, 3), c(-97.5)) 
  data3 <- calc_percentage(datos = ok_data[, "X3"])
  expect_equivalent(round(data3, 3), c(0)) 

  # file_data_one <- generate_mock_data_coor(i_year = 1901, ndata = 3) # 
  # ok_data <- calc_data_year(data = file_data_one$data)
  # ok_data[, "X1"] <- 1:length(ok_data[, "X1"])
  # ok_data[, "X2"] <- length(ok_data[, "X2"]):1
  # ok_data[, "X3"] <- 1

  # data1 <- calc_percentage(datos = ok_data[, "X1"])
  # expect_equivalent(round(data1, 3), c(14900)) 
  # data2 <- calc_percentage(datos = ok_data[, "X2"])
  # expect_equivalent(round(data2, 3), c(-99.333)) 
  # data3 <- calc_percentage(datos = ok_data[, "X3"])
  # expect_equivalent(round(data3, 3), c(0)) 

  file_data_one <- generate_mock_data_coor(i_year = 1871, ndata = 3) 
  ok_data <- calc_data_year(data = file_data_one$data)
  ok_data[, "X1"] <- 1:length(ok_data[, "X1"])
  ok_data[, "X2"] <- length(ok_data[, "X2"]):1
  ok_data[, "X3"] <- 1

  data1 <- calc_percentage(datos = ok_data[, "X1"])
  expect_equivalent(round(data1, 3), c(14900)) 
  data2 <- calc_percentage(datos = ok_data[, "X2"])
  expect_equivalent(round(data2, 3), c(-99.333)) 
  data3 <- calc_percentage(datos = ok_data[, "X3"])
  expect_equivalent(round(data3, 3), c(0)) 
})

#' testea la función coef_var
#'
#' @return None
#' @export
#'
test_that("coef_var", {
  file_data <- generate_mock_data_coor()
  ok_data <- calc_data_year(data = file_data$data)
  data <- coef_var(x = ok_data[, "X4"])
  expect_equivalent(round(data, 3), 0.016) # Función de Sergio, comprobamos que mantiene resultado
})

##########################################################

#' testea la función second_data_fill_data
#'
#' @return None
#' @export
#'
test_that("second_data_fill_data", {
  file_data <- generate_mock_data_coor()
  file_data$data[1:100, "X2"] <- NA
  file_data$data[140, "X3"] <- NA
  file_data$data[100, "X4"] <- NA
  file_data$data[200, "X5"] <- NA # Aug
  data <- second_data_fill_data(file_data = file_data)
  expect_equivalent(file_data$coor[c("X1", "X3", "X4", "X5"), ], data$coor)

  # Estación solitaria
  file_data$coor["X5", ] <- c(45, 90)
  data <- second_data_fill_data(file_data = file_data)
  expect_equivalent(file_data$coor[c("X1", "X3", "X4", "X5"), ], data$coor)
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
  expect_lt(data$data[1, 1], 0)

  # Estación solitaria
  file_data$data[12, "X5"] <- NA
  file_data$coor["X5", ] <- c(45, 90)
  data <- alexanderson_homogenize_data(file_data)
  expect_lt(data$data[1, 1], 0)
})

#' testea la función alexanderson_homogenize
#'
#' @return None
#' @export
#'
test_that("alexanderson_homogenize", {
  alexanderson_folder <- "test_result"

  file_data1 <- generate_mock_data_coor(ndata = 8)
  file_data1$data[c(1:100), 1] <- c(1:100)
  file_data1$data[c(101:dim(file_data1$data)[1]), 1] <- c(-101:-dim(file_data1$data)[1])
  data1 <- alexanderson_homogenize_data(file_data1)

  file_data2 <- generate_mock_data_coor(ndata = 8)
  file_data2$data[c(1:100), 1] <- c(101:200)
  file_data2$data[c(101:dim(file_data2$data)[1]), 1] <- c(-101:-dim(file_data2$data)[1])
  data2 <- alexanderson_homogenize_data(file_data2)

  data_result <- alexanderson_homogenize(data = list("file_data1" = file_data1, "file_data2" = file_data2), folder = alexanderson_folder)

  # Como todas las estaciones de data2 están en data1, sus valores finales tienen que ser los de data1 y no los de data2
  expect_equal(data_result$file_data1$data, data1$data)
  expect_equal(data_result$file_data2$data, data1$data)
})

#' testea la función calculate_statistics
#'
#' @return None
#' @export
#'
test_that("calculate_statistics", {
  file_data_one <- generate_mock_data_coor()
  file_data_one$data[is.na(file_data_one$data)] <- 35

  file_data_list <-  list("a1" = file_data_one, "a2" = file_data_one, "a3" = file_data_one)

  return_data <- calculate_statistics(data = file_data_list, data_ori = file_data_one$data)
  expect_equivalent(return_data$a1$regional_series[12], mean(file_data_one$data[12, ], na.rm = TRUE))

  # file_data_one <- generate_mock_data_coor(i_year = 1871, ndata = 3) #1901, 1931, 1951, 1961, 1981
  # file_data_one$data[, "X1"] <- 1:length(file_data_one$data[, "X1"])
  # file_data_one$data[, "X2"] <- length(file_data_one$data[, "X2"]):1
  # file_data_one$data[, "X3"] <- 1
  # return_data_1871 <- calculate_statistics_data(file_data = file_data_one, data_ori = file_data_one$data)
  # expect_equivalent(round(return_data_1961$magnitude_change_percentage$year, 4), c(7200, -98.6301, 0)) # MAL 35.13 -97.5 0

  # file_data_one <- generate_mock_data_coor(i_year = 1901, ndata = 3) #1901, 1931, 1951, 1961, 1981
  # file_data_one$data[, "X1"] <- 1:length(file_data_one$data[, "X1"])
  # file_data_one$data[, "X2"] <- length(file_data_one$data[, "X2"]):1
  # file_data_one$data[, "X3"] <- 1
  # return_data_1901 <- calculate_statistics_data(file_data = file_data_one, data_ori = file_data_one$data)
  # expect_equivalent(round(return_data_1961$magnitude_change_percentage$year, 4), c(7200, -98.6301, 0)) # MAL 35.13 -97.5 0

  # file_data_one <- generate_mock_data_coor(i_year = 1931, ndata = 3) #1901, 1931, 1951, 1961, 1981
  # file_data_one$data[, "X1"] <- 1:length(file_data_one$data[, "X1"])
  # file_data_one$data[, "X2"] <- length(file_data_one$data[, "X2"]):1
  # file_data_one$data[, "X3"] <- 1
  # return_data_1931 <- calculate_statistics_data(file_data = file_data_one, data_ori = file_data_one$data)
  # expect_equivalent(round(return_data_1961$magnitude_change_percentage$year, 4), c(7200, -98.6301, 0)) # MAL 35.13 -97.5 0
  

  # file_data_one <- generate_mock_data_coor(i_year = 1951, ndata = 3) #1901, 1931, 1951, 1961, 1981
  # file_data_one$data[, "X1"] <- 1:length(file_data_one$data[, "X1"])
  # file_data_one$data[, "X2"] <- length(file_data_one$data[, "X2"]):1
  # file_data_one$data[, "X3"] <- 1
  # return_data_1951 <- calculate_statistics_data(file_data = file_data_one, data_ori = file_data_one$data)
  # expect_equivalent(round(return_data_1961$magnitude_change_percentage$year, 4), c(7200, -98.6301, 0)) # MAL 35.13 -97.5 0

  # file_data_one <- generate_mock_data_coor(i_year = 1961, ndata = 3) #1901, 1931, 1951, 1961, 1981
  # file_data_one$data[, "X1"] <- 1:length(file_data_one$data[, "X1"])
  # file_data_one$data[, "X2"] <- length(file_data_one$data[, "X2"]):1
  # file_data_one$data[, "X3"] <- 1
  # return_data_1961 <- calculate_statistics_data(file_data = file_data_one, data_ori = file_data_one$data)
  # expect_equivalent(round(return_data_1961$magnitude_change_percentage$year, 4), c(7200, -98.6301, 0)) # MAL 35.13 -97.5 0 
  
  file_data_one <- generate_mock_data_coor(i_year = 1981, ndata = 3) #1901, 1931, 1951, 1961, 1981
  file_data_one$data[, "X1"] <- 1:length(file_data_one$data[, "X1"])
  file_data_one$data[, "X2"] <- length(file_data_one$data[, "X2"]):1
  file_data_one$data[, "X3"] <- 1
  return_data_1981 <- calculate_statistics_data(file_data = file_data_one, data_ori = file_data_one$data)
  expect_equivalent(round(return_data_1981$magnitude_change_percentage$year, 4), c(7200, -98.6301, 0)) # MAL 35.13 -97.5 0 
  # mobile_trends_data$year$X1$matriz_p # calculan mobile_trends_data c(1871, 1901, 1931)  
})


#' testea la función calculate_statistics_data
#'
#' @return None
#' @export
#'
test_that("calculate_statistics_data", {
  file_data <- generate_mock_data_coor()
  file_data$data[is.na(file_data$data)] <- 35
  data <- calculate_statistics_data(file_data = file_data, data_ori = file_data$data)
  expect_equivalent(data$regional_series[12], mean(file_data$data[12, ], na.rm = TRUE))
})

#' testea la función percentage_of_zeros
#'
#' @return None
#' @export
#'
test_that("percentage_of_zeros", {
  data <- c(1, 2, 3, 0, NA, NA, 0, 9, 0, 2, 0, 1)
  percentage <- percentage_of_zeros(data)
  expect_equivalent(percentage, 40)
})

#' testea la función delete_zero
#'
#' @return None
#' @export
#'
test_that("delete_zero", {
  file_data <- generate_mock_data_coor()

  file_data$data[c(12:16), "X1"] <- 0
  # percentage <- percentage_of_zeros(file_data$data[, "X1"])

  file_data$data[c(12:16), "X2"] <- 0
  file_data$data[c(20:dim(file_data$data)[1]), "X2"] <- 0
  # percentage <- percentage_of_zeros(file_data$data[, "X2"])

  data <- delete_zero(data = file_data$data)

  expect_equivalent(sum(!is.na(data[c(12:16), "X1"])), 0)
  expect_equivalent(sum(!is.na(data[c(12:16), "X2"])), 5)
})

#' testea la función save_delete_data
#'
#' @return None
#' @export
#'
test_that("save_delete_data", {
  folder_name <- "test_result"

  file_data <- generate_mock_data_coor()
  delete_data <- save_delete_data(file_data$data, file_data$data, folder_name)
  expect_equivalent(delete_data["X2", "data_delete"], 0)
  file_data_other <- file_data
  file_data_other$data[100, "X2"] <- NA
  file_data_other$data[102, "X2"] <- NA
  delete_data <- save_delete_data(file_data$data, file_data_other$data, folder_name)
  expect_equivalent(delete_data["X2", "data_delete"], 2)
})

#' testea la función delete_zones
#'
#' @return None
#' @export
#'
test_that("delete_zones", {
  file_data <- generate_mock_data_coor()
  file_data$coor[, "lat"] <- 40
  file_data$coor["X2", "lat"] <- 4

  data <- list("start" = file_data)
  delete_data <- delete_zones(data)

  expect_equivalent(rownames(delete_data$start$coor), c("X1", "X3", "X4", "X5"))
})

#' testea la función calculate_reconstruction_statistics
#'
#' @return None
#' @export
#'
test_that("calculate_reconstruction_statistics", {
  file_data <- generate_mock_data_coor()
  delete_data <- calculate_reconstruction_statistics(sim = file_data$data[, "X2"], obs = file_data$data[, "X2"])
  expect_equivalent(delete_data, c(1, 0, 0, 0))
})

#######################################################################################v

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

#   main_mediterranean_calculations(file_data = file_data, file_coor = file_coor, max_dist = max_dist_test)

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

#   main_mediterranean_calculations(file_data = file_data, file_coor = file_coor, max_dist = max_dist_test)

#   # sum(is.na(control_data$data[,"X5"]))

# })


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

