
library(R6)
library(lubridate)
library(dplyr)
library(mcga)

strata <- R6Class("strata",
                  public=list(
                    num_stratums=NULL,
                    num_partitions=NULL,
                    stratums=NULL,
                    partitions=NULL,
                    start_date=NULL,
                    end_date=NULL,
                    total_load=NULL,
                    rmse=NULL,
                    mse=NULL,
                    ci=NULL,
                    q_conc_df=NULL,
                    all_stratums_valid=NULL,
                    min_num_samples_per_stratum=3,
                    initialize=function( num_stratums, q_conc_df ) {
                      self$num_stratums <- num_stratums
                      self$num_partitions <- num_stratums - 1
                      self$start_date <- min( q_conc_df$date )
                      self$end_date <- max( q_conc_df$date )
                      self$q_conc_df <- q_conc_df
                      # 'strata' is a collection of 'stratums'
                      self$stratums <- vector("list", num_stratums )
                      self$partitions <- partitions$new( self$start_date, self$end_date, self$num_partitions )
                      for ( indx in 1:num_stratums ) {
                          self$stratums[[ indx ]] <- stratum$new()
                      }
                    },
                    all_valid=function() {
                      all_stratums_valid <- TRUE
                      for ( indx in 1:self$num_stratums ) {
                        if ( self$stratums[[ indx ]]$n_conc < self$min_num_samples_per_stratum ) {
                          all_stratums_valid <- FALSE
                        }
                      }
                      self$all_stratums_valid <- all_stratums_valid
                      return( all_stratums_valid )
                    },
                    calc_stats=function() {
                      self$rmse <- 0.0; self$mse <- 0.0; self$ci <- 0.0; self$total_load <- 0.0
                      for ( indx in 1:self$num_stratums ) {
                        self$total_load <- self$total_load + self$stratums[[ indx ]]$sl$load_total_corrected
                        self$mse        <- self$mse + self$stratums[[ indx ]]$sl$mse
                      }
                      if ( self$mse >= 0.0 ) {
                        self$rmse <- sqrt( self$mse )
                      } else {
                        self$rmse <- 1.0E+23
                      }
                    },
                    rearrange_stratums=function( partition_list ) {
                      self$rmse <- 0.0; self$mse <- 0.0; self$ci <- 0.0; self$total_load <- 0.0
                      num_partitions <- length(partition_list)
                      self$stratums[[1]]$update( self$start_date, partition_list[1], self$q_conc_df )
                      if ( num_partitions > 1 ) {
                        for ( i in 2:num_partitions ) {
                          self$stratums[[i]]$update(partition_list[i-1]+1, partition_list[i], self$q_conc_df  )
                        }
                      }
                      self$stratums[[self$num_stratums]]$update( partition_list[num_partitions] + 1, self$end_date, self$q_conc_df )
                    },
                    calc_loads=function() {
                      self$rmse <- 1E+23; self$mse <- 0.0; self$ci <- 1E+23; self$total_load <- 0.0
                      for ( indx in seq(1,self$num_stratums) ) {
                        self$stratums[[ indx ]]$calc_load()
                        self$total_load <- self$total_load + self$stratums[[ indx ]]$sl$load_total_corrected
                        self$mse        <- self$mse + self$stratums[[ indx ]]$sl$mse
                      }
                      if ( self$mse >= 0.0 & (!is.na( self$mse ) )) {
                        self$rmse <- sqrt( self$mse )
                      } else {
                        self$rmse <- 1.0E+23
                      }
                    },
                    update_rmse=function( x ) {
                      self$partitions$update( sort( x ) )
                      self$rearrange_stratums( self$partitions$get_date_values() )
                      if ( self$all_valid() ) {
                        self$calc_loads()
                      } else {
                        self$rmse <- 1E+23
                        self$mse  <- 1E+23
                      }
                      return( self$rmse )
                    }

                  ) )


