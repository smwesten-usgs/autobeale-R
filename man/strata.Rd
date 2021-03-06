% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strata.R
\docType{data}
\name{strata}
\alias{strata}
\title{Data structure to hold strata information.}
\format{An object of class \code{R6ClassGenerator} of length 24.}
\usage{
strata
}
\description{
In addition to storing information on the stratums and load, a number of
functions are associated with the data structure.
}
\section{Slots}{

\describe{
\item{\code{all_valid}}{returns \code{TRUE} if all underlying stratums contain sufficient concentration data}

\item{\code{calc_stats}}{loops over the underlying stratums to calculate stats for entire period}

\item{\code{rearrange_stratums}}{accepts a list of partition dates; moves partitions and recalculates stratum statistics}

\item{\code{calc_loads}}{loops over underlying stratums and calculates the load in each}

\item{\code{update_rmse}}{accepts new partitioning scheme, calculates the load, and returns the rmse}
}}

\examples{
# create new strata object with 3 stratums
number_of_stratums <- 3
mystrata <- strata$new( number_of_stratums, qcdata )

# create partitions to divide year into 3 roughly equal chunks (Jan-Apr, May-Aug, Sep-Dec)
mypartitions <- c( as.Date("1997-04-30"), as.Date("1997-08-31") )

# make call to rearrange stratum boundaries
# NOTE: this also copies the appropriate discharge and concentration data
#       into each stratum
mystrata$rearrange_stratums( mypartitions )

# now calculate loads for each stratum and for the strata overall
mystrata$calc_loads()
}
\keyword{datasets}
