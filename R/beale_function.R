#' Calculate pollutant loads by means of the Beale Ratio Estimator.
#'
#' This function is an R implementation of the code originally written
#' in Fortran. The Beale Ratio Estimator has been shown to provide an unbiased load estimate.
#'
#' Peter Richards, one of the authors of the Fortran implementation, describes it this way:
#'
#'    "The Beale Ratio Estimator is described in Tin (1965). It has been widely used for computation of
#'     substance loadings to receiving waters, particularly in the Great Lakes system. The initial code was
#'     developed by Ken Baun of Wisconsin DNR (1982). It was subsequently used and developed by Kevin McGunagle
#'     at the International Joint Commission (IJC), from whom the current version was received. It was then
#'     developed in the 1990’s by Pete Richards at the National Center for Water Quality Research
#'     at Heidelberg University (www.heidelberg.edu/academiclife/distinctive/ncwqr). Peter developed the code
#'     to work with the Macintosh computer of the time (“Motorola" architecture, which was replaced by
#'     Intel architecture ca. 2005), and included a Graphic User Interface (GUI) and an algorithm that objectively
#'     and sequentially searched out the best stratification given the data, with “best” being defined as
#'     minimizing the root-mean-square error (RMSE) given the data at hand. This is the current and final version.
#'     A Windows version of the code was also developed, with a minimal GUI."
#'
#' @param discharge_cms vector of daily discharge values given in cubic meters per second.
#' @param conc_mg_L paired vector of daily concentration values given in mg_L
#' @return A list with the following items:
#' \tabular{lll}{
#' Name \tab Description \cr
#' load_total_corrected \tab bias-corrected load for entire period of calculation \cr
#' load_daily_corrected \tab bias-corrected daily load \cr
#' load_mean \tab mean of daily loads over period of calculation \cr
#' load_daily \tab vector of individual daily loads \cr
#' confidence_interval \tab 95\% confidence interval about the total corrected load \cr
#' discharge_mean \tab mean discharge during calculation period \cr
#' discharge_mean_sample_days_only \tab mean discharge, including only days for which concentration data exist \cr
#' mse \tab mean squared error \cr
#' rmse \tab root mean square of errors \cr
#' df \tab number of degrees of freedom \cr
#' n_conc \tab number of concentration data points used in calculation \cr
#' n_discharge \tab number of discharge values used in calculation \cr
#' }
#' @examples
#' library(dplyr)
#' discharge_data <- read.table("ROCKFLOW.DAT", header=TRUE )
#' concentration_data <- read.table("ROCKNO3.DAT", header=TRUE )
#' qcdata <- left_join(discharge_data, concentration_data, by=c("date" = "date") )
#' head( qcdata )
#' beale( discharge_cmd=qcdata$discharge / 35.3146, conc_mg_L=qcdata$concentration )
#' @export
#' @md
beale <- function(discharge_cms, conc_mg_L) {

  discharge_mean         <- mean( discharge_cms, na.rm=TRUE )
  discharge_subset       <- discharge_cms[ !is.na(conc_mg_L) ]
  conc_subset            <- conc_mg_L[ !is.na(conc_mg_L) ]
  n_conc                 <- length( conc_subset )
  n_discharge            <- length(discharge_cms)
  load_daily             <- conc_subset * discharge_subset * 86.4
  load_mean              <- mean( load_daily,na.rm=TRUE )
  discharge_subset_mean  <- mean( discharge_subset, na.rm=TRUE )

  sql <- cov( discharge_subset, load_daily )     # covariance of flow and load: this is sxy in Pete's code; mass-vol/time
  sqq <- var( discharge_subset )                 # variance of discharge_cms: this is sx2 in Pete's code; (vol/time)^2
  sll <- var( load_daily )                       # variance of load: this is sy2 in Pete's code; mass^2

  numer                <- 1 + ( (1 / n_conc ) * sql / ( discharge_mean *load_mean ) )
  denom                <- 1 + ( ( 1 / n_conc ) * sqq / ( discharge_mean**2 ) )
  load_daily_corrected <- load_mean * ( discharge_mean /discharge_subset_mean ) * numer / denom
  load_total_corrected <- round( load_daily_corrected * n_discharge, 1 )

  ### Estimate the error
  sxy <- sql
  sx2 <- sqq
  sy2 <- sll
  sumx3 <- sum( ( discharge_subset - discharge_subset_mean )^3)
  sumx2y <- sum( ( discharge_subset - discharge_subset_mean )^2 * ( load_daily - load_mean ) )
  sumxy2 <- sum( ( discharge_subset - discharge_subset_mean ) * ( load_daily - load_mean )^2 )
  sx2y <- ( sumx2y / (n_conc - 1 ) ) / ( discharge_subset_mean^2 * load_mean )
  sxy2 <- ( sumxy2 / (n_conc - 1 ) ) / ( discharge_subset_mean * load_mean^2 )
  sx3  <- ( sumx3 / (n_conc - 1 ) ) / ( discharge_subset_mean^3 )

  s1 <- sx2 / ( discharge_subset_mean^2 )
  s2 <- sy2 / ( load_mean^2 )
  s3 <- sxy / ( discharge_subset_mean * load_mean )

  f1 <- ( 1 / n_conc ) * ( s1 + s2 - 2 * s3 )
  f2 <- ( 1 / n_conc )^2 * ( 2 * s1^2 - 4 * s1 * s3 + s3^2 + s1 * s2 )
  f3 <- ( ( 2 / n_conc ) / n_discharge ) * (sx3 - 2*sx2y + sxy2)

  mse <- (load_mean * discharge_mean / discharge_subset_mean )^2 * ( f1 + f2 + f3 ) * n_discharge^2
  rmse <- round((mse^0.5),3)

  df <- n_conc - 1

  # the following was in Rem's code:
  # CI <- qnorm(0.975)*(mse^0.5)/(df^0.5)

  CI <- ifelse( is.na( mse ), NA, qt( 0.975, df ) * sqrt( mse ) )

  returnval <- list(load_total_corrected=load_total_corrected,
                    load_daily_corrected=load_daily_corrected,
                    load_mean=load_mean,
                    load_daily=load_daily,
                    confidence_interval=CI,
                    discharge_mean=discharge_mean,
                    discharge_mean_sample_days_only=discharge_subset_mean,
                    mse=mse,
                    rmse=rmse,
                    df=df,
                    # sql=sql,
                    # sqq=sqq,
                    # sll=sll,
                    # f1=f1,
                    # f2=f2,
                    # f3=f3,
                    # s1=s1,
                    # s2=s2,
                    # s3=s3,
                    # sumx3=sumx3,
                    # sumx2y=sumx2y,
                    # sumxy2=sumxy2,
                    # sx2y=sx2y,
                    # sxy2=sxy2,
                    # sx3=sx3,
                    n_conc=n_conc,
                    n_discharge=n_discharge )
}
