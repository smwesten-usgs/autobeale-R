
library(R6)
library(lubridate)
library(dplyr)
library(mcga)

#' @title Data structure to hold partition information.
#'
#' @description In addition to storing information on the partitions, a number of
#' functions are associated with the data structure.
#' @slot initialize reserve space for the appropriate number of partitions
#' @slot get_date_fractions return the date fractions associated with the current set of partitions
#' @slot get_date_values return the dates associated with the current set of partitions
#' @export
#' @md
partitions <- R6Class("partitions",
                      public=list(
                        num_partitions=NULL,
                        start_date=NULL,
                        end_date=NULL,
                        pt=NULL,
                        initialize=function( start_date, end_date, num_partitions ) {
                          self$start_date <- start_date
                          self$end_date <- end_date
                          self$pt <- vector("list", num_partitions )
                          self$num_partitions <- num_partitions
                          for (i in 1:num_partitions) {
                            self$pt[[i]] <- partition$new()
                          }
                        },
                        update=function( date_fractions ) {
                          self$update( date_fractions )
                        },
                        update_fractions=function( date_fractions ) {
                          n <- length( date_fractions )
                          for (i in 1:n) {
                            pt <- self$pt[[i]]
                            pt$update( self$start_date, self$end_date, date_fractions[i] )
                          }
                        },
                        update_dates=function( date_values ) {
                          period_length <- as.numeric( self$end_date - self$start_date + 1 )
                          day_values <- as.numeric( date_values - self$start_date + 1 )
                          date_fractions <- day_values / period_length
                          self$update_fractions( date_fractions )
                        },
                        get_date_fractions=function() {
                          denominator <- as.numeric( self$end_date - self$start_date )
                          return_vals <- vector("numeric", self$num_partitions )
                          for( i in 1:self$num_partitions ) {
                            pt <- self$pt[[i]]
                            numerator <- as.numeric( pt$partition_date - self$start_date )
                            return_vals[i] <- numerator / denominator
                          }
                          return( return_vals )
                        },
                        get_date_values=function() {
                          return_vals <- structure( numeric( self$num_partitions ), class="Date" )
                          for( i in 1:self$num_partitions ) {
                            pt <- self$pt[[i]]
                            return_vals[i] <- pt$partition_date
                          }
                          return( return_vals )
                        }
                      ))

