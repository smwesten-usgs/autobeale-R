#' @title Convert cubic feet per second to cubic meters per second.
#'
#' @description Convert cubic feet per second to cubic meters per second.
#' @export
cfs_to_cms <- function( q_cfs ) {
  return( q_cfs / 35.3146665722 )
}
