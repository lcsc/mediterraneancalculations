% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions-mediterranean-calculations.r
\name{near_estations}
\alias{near_estations}
\title{Returns stations in order of proximity
Ignoring those that are more than 200 km away (NA in those cases)}
\usage{
near_estations(data, coor, max_dist)
}
\arguments{
\item{data}{monthly data}

\item{coor}{coordinates of the stations corresponding to the data}

\item{max_dist}{maximum distance between the series to be used}
}
\value{
correlation between stations
}
\description{
Returns stations in order of proximity
Ignoring those that are more than 200 km away (NA in those cases)
}
