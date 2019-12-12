#' Cel_devstages
#' 
#' Table of reference datasets and landmark developmental stages.
#' 
#' The timing of developmetnal stages were compiled from the WormBook \insertCite{wormbook}{wormRef} and the publications 
#' associated with the datasets used to build the references. 
#' This mainly concerns the larval reference \insertCite{hyun2013dampening}{wormRef}
#' 
#' @docType data
#' 
#' @section Format :
#'  A list with 2 dataframes:
#'  - `$devstages`
#'  - `$datasets`
#'  
#'  Both dataframes have the following fields:
#'  - `name`
#'  - `tstart`
#'  - `tend`
#'  - `tunit`
#'  
#'  
#' @references
#'    \insertAllCited{}
#' 
#' 
#' @importFrom Rdpack reprompt
"Cel_devstages"

