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
#' @param concentration_mg_L paired vector of daily concentration values given in mg_L
#' @return A list with the following items:
#' \tabular{lll}{
#' Name \tab Description \cr
#' load_total_corrected \tab bias-corrected load for entire period of calculation \cr
#' load_daily_corrected \tab bias-corrected daily load \cr
#' load_mean_daily \tab mean of daily loads over period of calculation \cr
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
#' flow_conc_df <- left_join(discharge_data, concentration_data, by=c("date" = "date") )
#' head( flow_conc_df )
#' beale( discharge_cms=cfs_to_cms( flow_conc_df$discharge ), concentration_mg_L=flow_conc_df$concentration )
#' @export
#' @md
beale <- function(discharge_cms, concentration_mg_L) {

  discharge_mean         <- mean( discharge_cms, na.rm=TRUE )
  discharge_subset       <- discharge_cms[ !is.na(concentration_mg_L) ]
  conc_subset            <- concentration_mg_L[ !is.na(concentration_mg_L) ]
  n_conc                 <- length( conc_subset )
  n_discharge            <- length(discharge_cms)
  load_daily             <- conc_subset * discharge_subset * 86.4
  load_mean_daily        <- mean( load_daily, na.rm=TRUE )
  discharge_subset_mean  <- mean( discharge_subset, na.rm=TRUE )

  theta <- 1 / ( n_conc - 1 )

  #sql <- cov( discharge_subset, load_daily )     # covariance of flow and load: this is sxy in Pete's code; mass-vol/time
  #sqq <- var( discharge_subset )                 # variance of discharge_cms: this is sx2 in Pete's code; (vol/time)^2
  #sll <- var( load_daily )                       # variance of load: this is sy2 in Pete's code; mass^2

  ### Estimate the error
  #sxy <- sql
  #sx2 <- sqq
  #sy2 <- sll

  #rTheta = 1_T_REAL / (REAL(n,kind=T_REAL)-1_T_REAL)
  #
  #do i=1,size(pConc%rConc)
  #
  #if(pConc(i)%iJulianDay >= pStratum%iStartDate &
  #   .and. pConc(i)%iJulianDay <= pStratum%iEndDate &
  #   .and. pConc(i)%lInclude) then
  #
  # update accumulators
  #
  # rSum_lq = rSum_lq + pConc(i)%rLoadTimesFlow &
  #   - (pStratum%rMeanSampleLoad * pStratum%rMeanSampleFlow)
  #
  # rSum_qq = rSum_qq + (pConc(i)%rFlow - pStratum%rMeanSampleFlow)**2
  #
  # rSum_ll = rSum_ll + (pConc(i)%rDailyLoad - pStratum%rMeanSampleLoad)**2
  #
  # rSum_q2l = rSum_q2l + (pConc(i)%rFlow - pStratum%rMeanSampleFlow)**2 &
  #   * (pConc(i)%rDailyLoad - pStratum%rMeanSampleLoad)
  #
  # rSum_q3 = rSum_q3 + (pConc(i)%rFlow - pStratum%rMeanSampleFlow)**3
  #
  # rSum_ql2 = rSum_ql2 + (pConc(i)%rFlow - pStratum%rMeanSampleFlow) &
  #   * (pConc(i)%rDailyLoad-pStratum%rMeanSampleLoad)**2
  #
  # end if
  #
  # end do
  #
  # pStratum%rS_lq = rTheta * rSum_lq     ! this is S_xy in Baun's paper, eqn B
  # pStratum%rS_qq = rTheta * rSum_qq     ! this is S_x**2 in Baun's paper, eqn B
  # pStratum%rS_ll = rTheta * rSum_ll
  # pStratum%rS_q2l = rTheta * rSum_q2l
  # pStratum%rS_q3 = rTheta * rSum_q3
  # pStratum%rS_ql2 = rTheta * rSum_ql2
  #

  sumSxy <- sum( ( discharge_subset - discharge_subset_mean ) * ( load_daily - load_mean_daily ) )
  sumSx2  <- sum( ( discharge_subset - discharge_subset_mean )^2 )
  sumSy2 <- sum( ( load_daily - load_mean_daily )^2 )
  sumSx3 <- sum( ( discharge_subset - discharge_subset_mean )^3)
  sumSx2y <- sum( ( discharge_subset - discharge_subset_mean )^2 * ( load_daily - load_mean_daily ) )
  sumSxy2 <- sum( ( discharge_subset - discharge_subset_mean ) * ( load_daily - load_mean_daily )^2 )

  Sxy <- n_conc * theta * sumSxy
  Sx2 <- n_conc * theta * sumSx2
  Sy2 <- n_conc * theta * sumSy2


  #rBiasCorr_numerator = 1_T_REAL + (1_T_REAL / REAL(pStratum%iNumSamples,kind=T_REAL) &
  #                                    * pStratum%rS_lq / rL_bar / rQ_bar)

  #rBiasCorr_denominator = 1_T_REAL + (1_T_REAL / REAL(pStratum%iNumSamples,kind=T_REAL) &
  #                                      * pStratum%rS_qq / (rQ_bar**2))


  numerator              <- 1 + theta * ( Sx2 / ( discharge_mean *load_mean_daily ) )
  denominator            <- 1 + theta *( Sy2 / ( discharge_mean**2 ) )
  bias_correction_factor <- numerator / denominator
  load_daily_biased      <- load_mean_daily * discharge_mean / discharge_subset_mean
  load_daily_corrected   <- load_daily_biased * bias_correction_factor
  load_total_biased      <- load_daily_biased * n_discharge
  load_total_corrected   <- load_daily_corrected * n_discharge
  bias_correction        <- load_total_corrected - load_total_biased

  Sx2y <-  theta * sumSx2y / ( discharge_subset_mean^2 * load_mean_daily )
  Sxy2 <-  theta * sumSxy2 / ( discharge_subset_mean * load_mean_daily^2 )
  Sx3  <-  theta * sumSx3  / ( discharge_subset_mean^3 )

  s1 <- Sx2 / ( discharge_subset_mean^2 )
  s2 <- Sy2 / ( load_mean_daily^2 )
  s3 <- Sxy / ( discharge_subset_mean * load_mean_daily )

  f1 <- theta * ( s1 + s2 - 2 * s3 )
  f2 <- theta^2 * ( 2 * s1^2 - 4 * s1 * s3 + s3^2 + s1 * s2 )
  f3 <- ( 2 * theta / n_discharge ) * (Sx3 - 2*Sx2y + Sxy2)

  #
  #
  # pStratum%rDailyBiasedLoadEstimate = pStratum%rMeanFlow * &
  #   pStratum%rMeanSampleLoad / pStratum%rMeanSampleFlow
  #
  # pStratum%rDailyCorrectedLoadEstimate = pStratum%rDailyBiasedLoadEstimate * &
  #   pStratum%rDailyLoadBiasCorrection
  #
  # rMSE_1 = rTheta * ( (pStratum%rS_qq / rQ_bar_sq) + &
  #                       (pStratum%rS_ll / rL_bar_sq) - &
  #                       (2_T_REAL * pStratum%rS_lq / (rL_bar * rQ_bar)) )
  #
  # rMSE_2 = rTheta**2 *(  2_T_REAL*(pStratum%rS_qq**2/rQ_bar_sq**2) - &
  #                          (4_T_REAL*pStratum%rS_qq*pStratum%rS_lq/(rQ_bar_sq * rL_bar * rQ_bar)) + &
  #                          (pStratum%rS_lq**2/((rL_bar*rQ_bar)**2)) + &
  #                          ((pStratum%rS_qq * pStratum%rS_ll)/(rQ_bar_sq * rL_bar_sq))  )
  #
  #
  # ! this term is from Tin (1965), equation V(t2), top of pp 299.
  # rMSE_3 = ( 2_T_REAL * rTheta /   REAL(pStratum%iNumDays,kind=T_REAL) ) * &
  #   ( (pStratum%rS_q3 / rQ_bar**3) - &
  #       (2_T_REAL*pStratum%rS_q2l / (rL_bar*rQ_bar**2)) + &
  #       (pStratum%rS_ql2 / (rL_bar**2 *rQ_bar)) )
  #
  #







  mse <- (load_mean_daily * discharge_mean / discharge_subset_mean )^2 * ( f1 + f2 + f3 ) * n_discharge^2
  rmse <- round((mse^0.5),3)

  df <- n_conc - 1

  # the following was in Rem's code:
  # CI <- qnorm(0.975)*(mse^0.5)/(df^0.5)

  CI <- ifelse( is.na( mse ), NA, qt( 0.975, df ) * sqrt( mse ) )

  returnval <- list(load_total_corrected=load_total_corrected,
                    load_daily_corrected=load_daily_corrected,
                    load_mean_daily=load_mean_daily,
                    load_daily=load_daily,
                    load_daily_biased=load_daily_biased,
                    load_daily_corrected=load_daily_corrected,
                    load_total_biased=load_total_biased,
                    load_total_corrected=load_total_corrected,
                    bias_correction=bias_correction,
                    confidence_interval=CI,
                    discharge_mean=discharge_mean,
                    discharge_mean_sample_days_only=discharge_subset_mean,
                    mse=mse,
                    rmse=rmse,
                    df=df,
                    conc_subset=conc_subset,
                    discharge_subset=discharge_subset,
                    Sxy=Sxy,
                    Sy2=Sy2,
                    Sx2=Sx2,
                    f1=f1,
                    f2=f2,
                    f3=f3,
                    s1=s1,
                    s2=s2,
                    s3=s3,
                    sumSx3=sumSx3,
                    sumSx2y=sumSx2y,
                    sumSxy2=sumSxy2,
                    Sx2y=Sx2y,
                    Sxy2=Sxy2,
                    Sx3=Sx3,
                    n_conc=n_conc,
                    n_discharge=n_discharge )
}
