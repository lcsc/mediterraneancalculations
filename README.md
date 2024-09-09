To install the package, it can be done directly in R (we recommend upgrading to the latest available version beforehand) with the instructions:
```r
install.packages("devtools")
library(devtools)
install_local('mediterraneancalculations_0.2-2.tar.gz') # This file has to be in the R directory
```
More installation ways can be seen here https://riptutorial.com/r/example/5556/install-package-from-local-source


The package code has been posted on https://github.com/lcsc/mediterraneancalculations And it can be installed from Github as follows
```r
install.packages("devtools")
library(devtools)
install_github("lcsc/mediterraneancalculations")
```


Once the package is installed, please load the following:
```r
library(mediterraneancalculations)
```

We call the main function:
```r
main_mediterranean_calculations(file_data, file_coor)
```

Where the parameters are:
1. file_data: CSV file of station data. # This file has to be in the R directory
2. file_coor: CSV file of station coordinates. # This file has to be in the R directory
In the case of the example files, the usage would be:
```r
result <- main_mediterranean_calculations("israel_data.csv", "israel_coor.csv")
```

It is very important that the input data is formatted like the sample data as follows:
* Coordinate file, rows separated by ";" and using "." to indicate decimals. It has to have 3 columns, "coor", "lat" and "lon". In "coor" we will have the id of the stations (if they are numbers, better with an X in front) that must coincide for the ids that we use in the data file. In "lat" and "lon" we will have the coordinates in latitude, longitude in WGS84 https://spatialreference.org/ref/epsg/wgs-84/ Please you follow exactly the format of the example files of the Israel dataset.
* Data file, with columns separated from rows by ";" and using "." to indicate decimals. The first columns are "year" and "month" and then the station IDs will come. So each row will have a year, followed by a month (identified by integers from 1 to 12) and then all the station data will be available. Missing data is expressed as NA or a blank.

In order to optimise the long-term series that generate the software, we recommend you to generate manually some long-term series (if possible according to your data availability). In practical terms, if you have some precipitation stations that cover different periods and that are located at short distances (usually in the same city), we would suggest you merge them in one single series before including the series in the database. For this purpose, you should use as the coordinates of the merged series those of the last precipitation station with data available. 

This should be done exclusively for the long-term precipitation series. Please, you do not expend time to reconstruct manually the precipitation series in all of the cases (including short series) since the software will do it. This initial manual approach should be exclusively done if in your country some long-term data is available.

This is what we did with the data of some long-term stations in Spain. As a representative example, in Zaragoza (our city), the first meteorological station was installed in 1858 and it collected data to 1905, when it moved to other location in the city. This station collected data to 1945, when it moved to the current location. Distances were short between observatories (< 5km.) so we directly generated a single series that started in 1858 and finished in 2020, merging the series in the three available observatories. Otherwise, if we would consider these three series independently, the software would not reconstruct any long-term series as they do not contain the sufficient required original data in each station to generate a continuous long term series.

Please, you try to follow this manual approach in the case of the series that start before 1930 to check if you can generate more stations that contain long-term records (maybe considering maximum distances of 10-15 km between series). We would not expend that tisk task requires a lot of time of effort for you as not many stations are expected to be available data for the long term.  

For example, in Spain we have generated 28 complete stations from 1871 and 93 from 1901 thanks to follow this approach. 

The results should be, in the "results" folder that will contain 5 data files:
* data_1871.csv: Filled and homogenised data of the stations that have been completed from the year 1871 to 2020.
* data_1901.csv: Filled and homogenised data of the stations that have been completed from the year 1901 to 2020.
* data_1931.csv: Filled and homogenised data of the stations that have been completed from the year 1931 to 2020.
* data_1951.csv: Filled and homogenised data of the stations that have been completed from the year 1951 to 2020.
* data_1961.csv: Filled and homogenised data of the stations that have been completed from the year 1961 to 2020.
* data_1981.csv: Filled and homogenised data of the stations that have been completed from the year 1981 to 2020.
And it will also contain a coordinate file for each data file.

A mediterranean_calculations.rds file will also be generated with the necessary data for the study. This file can be loaded with:
```r
data_save <- readRDS("mediterranean_calculations.rds")
```

The content of the files can be explored for each of the result groups (using seasonal and annual data 1871, 1901, 1931, 1951, 1961 and 1981):
* coor: Station Coordinates
* data_change: The number of data that have been removed in the quality control analysis. 
* regional_series: Average country precipitation series
* pvalue_trend: The p-value of the trend analysis
* pvalue_corrected_trend: The p-value corrected according to the method discussed in Montpellier.
* magnitude_change_mm: Magnitude of change in precipitation (in mm.)
* magnitude_change_percentage: Magnitude of change in precipitation (in %)
* average_precipitation_stations:  The mean precipitation in each station (annual seasonal) 
* cv_precipitation_stations: The coefficient of variation of precipitation in each station (annual and seasonal)
* change_values_homog: The coefficients used in homogenization in each station.
* break_points_homog: The break points of inhomogeneities
* pvalue_trend_homog: p-values of the series before homogeneity testing and correction.
* magnitude_change_mm_homog: Magnitude of change in mm in the series before homogeneity testing and correction.
* magnitude_change_percentage_homog: Magnitude of change in % in the series before homogeneity testing and correction.
* statistics_data_removed: Description of the origin of the data in the final series (original or reconstructed).
* availability_data_original: Coordinates with original (0 or 1) data in each month.

A sample code to open these RDS files, for example, the ones from [Israel](https://zenodo.org/records/10022618/files/Israel.rds): 
```r
results_israel <- readRDS("Israel.rds")

# Identify the different results groups
names(results_israel)

# Select the different groups, for example the coordinates of the stations with data since 1981.
coordinates <- results_israel$start_1871$coor

# Magnitude of change in mm in spring for stations with data since 1951.
Spring_change <- results_israel$start_1951$magnitude_change_mm$spring
```
