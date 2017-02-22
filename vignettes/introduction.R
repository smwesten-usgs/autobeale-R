## ---- initialization, warning=FALSE,message=FALSE------------------------

library(autobeale)
library(dplyr)
library(mcga)



## ---- results='asis'-----------------------------------------------------

discharge_file     <- system.file( "extdata", "ROCKFLOW.DAT", package="autobeale")
concentration_file <- system.file( "extdata", "ROCKNO23.DAT", package="autobeale")

discharge_data     <- read.table( discharge_file, header=TRUE )
concentration_data <- read.table( concentration_file, header=TRUE )
qcdata <- dplyr::left_join(discharge_data, concentration_data, by=c("date" = "date") )

# don't forget to convert the discharge value to cubic meters per second!
qcdata$discharge <- cfs_to_cms( qcdata$discharge )


## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(qcdata, 10))

## ---- unstratified_beale, results='asis'---------------------------------
result_list <- beale( discharge_cms=qcdata$discharge / 35.3146, conc_mg_L=qcdata$concentration )
print( result_list )

## ---- date munging-------------------------------------------------------

library(lubridate)

qcdata$date <- lubridate::ymd( qcdata$date )


## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(qcdata, 10))

## ---- strata and partitioning--------------------------------------------
number_of_stratums <- 5
mystrata <- strata$new( number_of_stratums, qcdata )

## ---- genetic_algorithm--------------------------------------------------

 m <- mcga( popsize=200,
            chsize=4,
            minval=0.0,
            maxval=1.0,
            elitism=10,
            maxiter=40,
            crossprob=1.0,
            mutateprob=0.05,
            evalFunc=mystrata$update_rmse )

 cat("Best chromosome:\n")
 cat(m$population[1,],"\n")
 cat("Cost: ",m$costs[1],"\n")

# create a copy of the strata object
best_strata <- mystrata$clone()

# rerun with the optimized parameter
best_strata$partitions$update( sort( m$population[ 1, ]) )
best_strata$rearrange_stratums( best_strata$partitions$get_date_values() )

best_strata$calc_loads()
results_df <- best_strata$summarize_results()

# round values to produce a more readable table
results_df <- data.frame(lapply(results_df, function(y) if(is.numeric(y)) round(y, 2) else y)) 

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable( results_df, digits=2 )

