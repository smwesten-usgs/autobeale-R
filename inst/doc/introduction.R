## ---- initialization, warning=FALSE,message=FALSE------------------------

library(autobeale)
library(dplyr)



## ---- results='asis'-----------------------------------------------------

discharge_file     <- system.file( "extdata", "ROCKFLOW.DAT", package="autobeale")
concentration_file <- system.file( "extdata", "ROCKNO23.DAT", package="autobeale")

discharge_data     <- read.table( discharge_file, header=TRUE )
concentration_data <- read.table( concentration_file, header=TRUE )
qcdata <- dplyr::left_join(discharge_data, concentration_data, by=c("date" = "date") )


## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(qcdata, 10))

## ---- unstratified_beale, results='asis'---------------------------------
result_list <- beale( discharge_cms=qcdata$discharge / 35.3146, conc_mg_L=qcdata$concentration )
print( result_list )

