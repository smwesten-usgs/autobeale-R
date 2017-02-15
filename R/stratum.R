
library(R6)
library(lubridate)
library(dplyr)
library(mcga)


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
                       idx <- qcdata$date >= start_date & qcdata$date <= end_date
                       self$q <- qcdata$q[ idx ]
                       self$c <- qcdata$c[ idx ]
                       self$n_conc <- sum( !is.na( self$c ) )
                     },
                     calc_load=function() {
                       self$sl <- beale(discharge_cms=self$q, conc_mg_L=self$c)

#                       self$q_mean <- sl$q_mean
#                       self$c_mean <- sl$c_mean
#                       self$l_mean <- sl$l_mean
#                       self$n_discharge <- sl$n_discharge
#                       self$n_conc <- sl$n_conc
                     }

                   ))

