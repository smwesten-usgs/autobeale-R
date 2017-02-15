
library(R6)
library(lubridate)
library(dplyr)
library(mcga)


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
                          n <- length( date_fractions )
                          for (i in 1:n) {
                            pt <- self$pt[[i]]
                            pt$update( self$start_date, self$end_date, date_fractions[i] )
                          }
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

