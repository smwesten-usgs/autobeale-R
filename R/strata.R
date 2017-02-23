
library(R6)
library(lubridate)
library(dplyr)
library(mcga)

#' @title Data structure to hold strata information.
#'
#' @description In addition to storing information on the stratums and load, a number of
#' functions are associated with the data structure.
#' @slot all_valid returns `TRUE` if all underlying stratums contain sufficient concentration data
#' @slot calc_stats loops over the underlying stratums to calculate stats for entire period
#' @slot rearrange_stratums accepts a list of partition dates; moves partitions and recalculates stratum statistics
#' @slot calc_loads loops over underlying stratums and calculates the load in each
#' @slot update_rmse accepts new partitioning scheme, calculates the load, and returns the rmse
#' @examples
#' # create new strata object with 3 stratums
#' number_of_stratums <- 3
#' mystrata <- strata$new( number_of_stratums, qcdata )
#'
#' # create partitions to divide year into 3 roughly equal chunks (Jan-Apr, May-Aug, Sep-Dec)
#' mypartitions <- c( as.Date("1997-04-30"), as.Date("1997-08-31") )
#'
#' # make call to rearrange stratum boundaries
#' # NOTE: this also copies the appropriate discharge and concentration data
#' #       into each stratum
#' mystrata$rearrange_stratums( mypartitions )
#'
#' # now calculate loads for each stratum and for the strata overall
#' mystrata$calc_loads()
#' @export
#' @md
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
                    df=NULL,
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
                      self$rmse <- 0.0; self$mse <- 0.0; self$ci <- 0.0; self$df <- 1E+23;  self$total_load <- 0.0
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
                      self$rmse <- 0.0; self$mse <- 0.0; self$ci <- 0.0; self$df <- 1E+23;  self$total_load <- 0.0
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
                      self$rmse <- 1E+23; self$mse <- 0.0; self$ci <- 1E+23; self$df <- 1E+23; self$total_load <- 0.0
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
                    },
                    summarize_results=function() {
                      n <- self$num_stratums + 1
                      df <- data.frame( stratum=character(n), from_date=structure(numeric(n), class="Date"),
                                        to_date=structure( numeric(n), class="Date"),
                                        n_conc=numeric(n), n_q=numeric(n), q_mean=numeric(n),
                                        q_sample_mean=numeric(n), load_daily_biased=numeric(n),
                                        load_daily_corrected=numeric(n), mse=numeric(n),
                                        rmse=numeric(n), df=numeric(n), load_total=numeric(n), ci=numeric(n),
                                        bias_correction=numeric(n),
                                        stringsAsFactors=FALSE )
                      for ( i in 1:self$num_stratums ) {
                        sl <- self$stratums[[i]]$sl
                        df$stratum[i] <- as.character( x=i )
                        df$from_date[i] <- self$stratums[[i]]$start_date
                        df$to_date[i] <- self$stratums[[i]]$end_date
                        df$n_conc[i] <- sl$n_conc
                        df$n_q[i] <- sl$n_discharge
                        df$q_mean[i] <- sl$discharge_mean
                        df$q_sample_mean[i] <- sl$discharge_mean_sample_days_only
                        df$load_daily_biased[i] <- sl$load_daily_biased
                        df$load_daily_corrected[i] <- sl$load_daily_corrected
                        df$bias_correction[i] <- sl$bias_correction
                        df$mse[i] <- sl$mse
                        df$rmse[i] <- sl$rmse
                        df$df[i] <- sl$df
                        df$load_total[i] <- sl$load_total_corrected
                        df$ci[i] <- sl$confidence_interval
                      }
                      df$stratum[n] <- "all strata"
                      df$from_date[n] <- self$stratums[[1]]$start_date
                      df$to_date[n] <- self$stratums[[i]]$end_date
                      df$n_conc[n] <- sum( !is.na( self$q_conc_df$concentration ) )
                      df$n_q[n] <- sum( !is.na( self$q_conc_df$discharge ) )
                      df$q_mean[n] <- mean( self$q_conc_df$discharge, na.rm=TRUE )
                      df$q_sample_mean[n] <- mean( self$q_conc_df$discharge[ !is.na( self$q_conc_df$concentration ) ], na.rm=TRUE )
                      df$load_daily_biased[n] <- mean( df$load_daily_biased[1:n-1] )
                      df$load_daily_corrected[n] <- mean( df$load_daily_corrected[1:n-1] )
                      df$bias_correction[n] <- df$load_daily_corrected[n] - df$load_daily_biased[n]
                      df$mse[n] <- self$mse
                      df$rmse[n] <- self$rmse
                      df$df[n] <- self$df
                      df$load_total[n] <- sum( df$load_total[1:n-1])
                      df$ci[n] <- self$ci

                      return( df )
                    }

                ) )
