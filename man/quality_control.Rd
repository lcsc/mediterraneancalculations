% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions-mediterranean-calculations.r
\name{quality_control}
\alias{quality_control}
\title{Control de calidad
Estaciones con menos de 20 años de datos retirar
Usando las 10 más cercanas a menos de 200 km, desechar si promedio de percentil se diferencia en más de 0.6 
o en más de 0.5 para datos 0}
\usage{
quality_control(data, coor, max_dist)
}
\arguments{
\item{data}{datos}

\item{coor}{coordenadas}

\item{max_dist}{máxima distancia entre 2 estaciones para usar una para evaluar o completar la otra}
}
\value{
data y coor con los datos que no pasan el control eliminados
}
\description{
Control de calidad
Estaciones con menos de 20 años de datos retirar
Usando las 10 más cercanas a menos de 200 km, desechar si promedio de percentil se diferencia en más de 0.6 
o en más de 0.5 para datos 0
}
