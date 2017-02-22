
library(R6)
library(lubridate)
library(dplyr)
library(mcga)


#' @title Data structure to hold information for a single stratum.
#'
#' @description Data structure to holds load and concentration data for
#' an individual stratum, and also makes the call to `beale` to estimate
#' the stratum's load.
#' @slot update adjusts the bounds of the current stratum based on the input
#' start and end date; subsets appropriate discharge and concentration data
#' @slot calc_load calculate the load for the stratum defined by this data structure
#' @export
#' @md
stratum <- R6Class("stratum",
                   public=list(
                     start_date=NULL,
                     end_date=NULL,
                     n_conc=NULL,
                     sl=NULL,
                     q=NULL,
                     c=NULL,
                     update=function(start_date, end_date, qcdata ) {
                       self$start_date <- start_date
                       self$end_date <- end_date
                       indx <- qcdata$date >= start_date & qcdata$date <= end_date
                       self$q <- qcdata$discharge[indx]
                       self$c <- qcdata$concentration[indx]
                       self$n_conc <- sum( !is.na( self$c ) )
                     },
                     calc_load=function() {
                       self$sl <- beale(discharge_cms=self$q, conc_mg_L=self$c)
                     }

                   ))

