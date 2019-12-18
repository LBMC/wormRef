#' Cel_larval
#' 
#' Transcriptomic data from *C. elegans* larval samples.
#' 
#' Three time series are included in this dataset.
#' 
#'  - `O.20`: \insertCite{hyun2013dampening}{wormRef}
#'  - `O.25`: \insertCite{hyun2013dampening}{wormRef}
#'  - `H`: \insertCite{hendriks2014extensive}{wormRef}
#' 
#' 
#' The gene expression matrix is \eqn{log(X + 1)} of quantile-normalized raw data.
#' 
#' @docType data
#' 
#' @eval data_format()
#' 
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49043}{GSE49043} 
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52861}{GSE52861}
#' 
#' @references
#'    \insertAllCited{}
#' 
#' 
#' @importFrom Rdpack reprompt
"Cel_larval"
