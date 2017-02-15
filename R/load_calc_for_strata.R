
library(R6)
library(mcga)


load_calc_for_strata <- function( strata_object ) {

  m <- mcga( popsize=300,
             chsize=strata_object$num_stratums - 1,
             minval=0.0,
             maxval=1.0,
             elitism=10,
             maxiter=60,
             crossprob=1.0,
             mutateprob=0.05,
             evalFunc=strata_object$update_rmse )

  cat("Best chromosome:\n")
  print(m$population[1,])
  cat("Cost: ",m$costs[1],"\n\n")

  best_strata <- strata_object$clone()

  # rerun with the optimized parameter
  best_strata$partitions$update( sort( m$population[ 1, ]) )
  best_strata$rearrange_stratums( best_strata$partitions$get_date_values() )

  best_strata$calc_loads()

  print( best_strata )
  print( best_strata$stratums )

}

