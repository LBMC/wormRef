#' Cel_YA_1
#' 
#' Transcriptomic data from *C. elegans* aging adults/ rrf-3 mutants .
#' 
#' One (unpublished) time series is included in this dataset.
#' 
#'  - Byrne et al. 
#' 
#' 
#' The gene expression matrix is \eqn{log(X + 1)} of quantile-normalized raw data, 
#' filtered to keep genes with a strong aging signal.
#' 
#' @docType data
#' 
#' @eval data_format()
#' 
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93826}{GSE93826} 
#' 
#' @references
#'    \insertAllCited{}
#' 
#' 
#' @importFrom Rdpack reprompt
"Cel_aging_1"
