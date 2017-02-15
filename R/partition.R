
library(R6)
library(lubridate)
library(dplyr)
library(mcga)


partition <- R6Class("partition",
                     public=list(
                       partition_date=NULL,
                       date_fraction=NULL,
                       update=function( start_date, end_date, date_fraction ) {
                         self$date_fraction <- date_fraction
                         deltat <- end_date - start_date
                         self$partition_date <- start_date + deltat * date_fraction
                       }

                     ))

