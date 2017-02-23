library(R6)
library(lubridate)
library(dplyr)
library(mcga)


#' @title Data structure to hold information for a Beale calculation.
#'
#' @description Data structure to store the variance terms and bias correction
#' factors that go into the Beale ratio estimator.
#'
#' Definition of variables used:
#' -----------------------------
#'
#' L     :: mean pollutant load, L (kg), or loading rate Ld, (kg/day)
#' Q     :: mean flow (cubic meters per day)
#' q     :: mean flow at times of sampling (cubic meters per day)
#' l     :: mean flow-weighted loading rate (kg/day)
#' B     :: bias correction factor
#' n     :: number of samples
#' x     :: discrete flow measurement, (l/sec)
#' c     :: discrete concentration measurement (mg/l)
#' y     :: discrete loading estimate, x * c (mg/sec)
#' MSE   :: mean-square-error
#' RMSE  :: square root of mean-square-error
#' N     :: duration (days)
#' edf   :: effective degrees of freedom over strata
#' CI    :: confidence interval
#'
#'
#'
#'      if(r(i) .gt. 1.0d0) then
#' sxy=(sumfl(i)-(r(i)*avflow*avload))/(r(i)-1.0d0)
#' sx2=(sumf2(i)-(r(i)*avflow*avflow))/(r(i)-1.0d0)
#' sy2=(suml2(i)-(r(i)*avload*avload))/(r(i)-1.0d0)
#' if (avload.gt.0.0d0) then
#' sx2y=(sumx2y(i)/(r(i)-1.0d0))/(avflow**2*avload)
#' sxy2=(sumxy2(i)/(r(i)-1.0d0))/(avflow*avload**2)
#' sx3= (sumx3(i)/(r(i)-1.0d0))/(avflow**3)
#'
#' SUBSCRIPTS:
#' -----------
#' h     :: stratum identifier
#' d     :: daily estimate
#' @slot calc_load calculate the load for the stratum defined by this data structure
#' @export
#' @md
beale <- R6Class("beale",
                   public=list(
                     q=NULL,
                     c=NULL,
                     date=NULL,
                     c_date=NULL,
                     initialize=function( discharge_cms, concentration_mg_l ) {
                       self$c                       <- concentration_mg_l[ !is.na( concentration_mg_l ) ]
                       self$discharge_X             <- discharge_cms
                       self$sample_discharge_x_i    <- discharge_cms[ !is.na( concentration_mg_l ) ]
                       self$num_samples_n_h         <- length( self$c )
                       self$duration_N              <- length( self$discharge_X )
                       self$mean_sample_discharge_q <- mean( self$sample_discharge_x_i, na.rm=TRUE )
                       self$mean_discharge_Q        <- mean( self$discharge_X, na.rm=TRUE )
                       self$theta                   <- 1 / ( self$num_samples_n_h )
                       self$theta_n_minus_1         <- 1 / ( self$num_samples_n_h - 1 )
                     }
                     calculate_discrete_loading_estimates=function() {
                       self$discrete_load_y_i <- self$sample_discharge_x_i * self$c * 86.4
                     },
                     calculate_mean_flow_weighted_loading_rate=function() {
                       self$l_h <- self$theta_n_minus_1 * sum( self$discrete_load_y_i )
                     },
                     calculate_bias_factor_numerator=function() {

                     },
                     calculate_bias_factor_denominator=function() {

                     },
                     #
                     # do (l=1,ndays)
                     #   if (stratum(l).ne.i) cycle
                     #   nf(i)=nf(i)+1
                     #   flowmu(i)=flowmu(i)+flow(l)
                     #     if (conc(l) .ne. -1.0d0) then
                     #       templ=flow(l)*conc(l)*86.4d0
                     #       r(i)=r(i)+1.0d0
                     #       suml(i)=suml(i)+templ
                     #       suml2(i)=suml2(i)+templ*templ
                     #       sumf(i)=sumf(i)+flow(l)
                     #       sumf2(i)=sumf2(i)+flow(l)*flow(l)
                     #       sumfl(i)=sumfl(i)+flow(l)*templ
                     #       fl1=flow(l)/0.028316849d0                  !put it also in cfs for listing
                     #       if (bigio) write (2,3) dates(l),templ,flow(l),fl1,conc(l)/(concfac*loadfac)
                     #     end if
                     #   repeat
                     #   if (r(i).le.1.0d0) then      !insufficient data to calculate a load for stratum
                     #     if (bigio) write (2,6)  oldtribname,strat(i,5)
                     #     write (3,6) oldtribname,strat(i,5)
                     #     needs_work=.true.
                     #     go to 300
                     #   end if
                     #   if (nf(i).le.0) then            !no flow data - can't calculate a load
                     #     if (bigio) write (2,5)
                     #     write (3,5)
                     #     call showalert(3, 'No flow data in one or more strata!'//
                     #    &          'Please adjust the strata...')
                     #     go to 300
                     #   end if
                     #   fpc=1.0d0-(r(i)/dble(nf(i)))
                     #   fpc=1.0d0
                     #   flowmu(i)=flowmu(i)/dble(nf(i))
                     #   avflow=sumf(i)/r(i)
                     #   avload=suml(i)/r(i)
                     #
                     #
                     # do (l=1,ndays)        !now calculate the third order terms
                     #   if (stratum(l).ne.i) cycle
                     #   if (conc(l).ne.-1.0d0) then
                     #     sumx3(i)=sumx3(i) + (flow(l)-avflow)**3
                     #     sumx2y(i)=sumx2y(i) + (flow(l)-avflow)**2 * (flow(l)*conc(l)*86.4d0-avload)
                     #     sumxy2(i)=sumxy2(i) + (flow(l)-avflow) * (flow(l)*conc(l)*86.4d0-avload)**2
                     #   end if
                     # repeat
                     #
                     calculate_moments=function() {
                       qh_lh   <- self$mean_sample_discharge_q * self$l_h
                       lh2     <- self$l_h**2
                       qh2     <- self$mean_sample_discharge_q**2
                       qh4     <- qh2 * qh2
                       qh2_lh2 <- qh2 * lh2
                       #' subpart of equation B, Baun (1982); Baun drops the "nh - 1" term in the publication
                       # Peter's code: sxy=(sumfl(i)-(r(i)*avflow*avload))/(r(i)-1.0d0)
                       self$Sxy <- self$theta_n_minus_1 * ( sum( self$sample_discharge_x_i * self$discrete_load_y_i )
                                       - self$num_samples_n_h * qh_lh )
                       #' subpart of equation B, Baun (1982); Baun appears to mean X_hi rather than Y_hi
                       #' in definition for (Sx_h)^2
                       #' Peter's code: sx2=(sumf2(i)-(r(i)*avflow*avflow))/(r(i)-1.0d0)
                       self$Sx2 <- self$theta_n_minus_1 * ( sum( self$sample_discharge_x_i**2 )
                                       - self$num_samples_n_h * qh2 )
                       #' Peter's code: sy2=(suml2(i)-(r(i)*avload*avload))/(r(i)-1.0d0)
                       self$Sy2 <- self$theta_n_minus_1 * ( sum( self$sample_discharge_x_i**2 )
                                       - self$num_samples_n_h * lh2 )
                       #self$sumx3 <- sum( ( self$sample_discharge_x_i - self$mean_sample_discharge_q )**3 )
                       numerator   <- 1 + self$theta * self$Sxy / qh_lh
                       denominator <- 1 + self$theta * self$Sx2 / qh2
                       #' equation B, Baun (1982)
                       self$bias_correction_factor_B <- numerator / denominator

                     }

                   ))
