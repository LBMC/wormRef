#' Cel_embryo
#' 
#' Transcriptomic data from *C. elegans* embryonic samples.
#' 
#' One time series is included in this dataset.
#'  - `L`: \insertCite{levin2016mid}{wormRef}
#' 
#' 
#' The gene expression matrix is \eqn{log(X + 1)} of quantile-normalized raw data.
#' 
#' @docType data
#' 
#' @eval data_format()
#' 
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60755}{GSE60755}
#' 
#' @references
#'    \insertAllCited{}
#' 
#' 
#' @importFrom Rdpack reprompt
"Cel_embryo"
