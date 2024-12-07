% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediterranean-calculations.r
\name{second_data_fill_data}
\alias{second_data_fill_data}
\title{Performs a second fill
For each dataset (1870, 1900...)
Stations with more than 90 or 95% of filled data
Sort the stations by correlation (minimum 0.5) and fill using the 10 methods...
Stations without total fill are discarded}
\usage{
second_data_fill_data(file_data, fillable_years = 36, max_dist = NA)
}
\arguments{
\item{file_data}{path to the data file}

\item{fillable_years}{years that can be filled with the station's monthly average}

\item{max_dist}{maximum allowed distance for filling}
}
\value{
None
}
\description{
Performs a second fill
For each dataset (1870, 1900...)
Stations with more than 90 or 95% of filled data
Sort the stations by correlation (minimum 0.5) and fill using the 10 methods...
Stations without total fill are discarded
}
