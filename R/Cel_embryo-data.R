#' Cel_embryo
#' 
#' Transcriptomic data from *C. elegans* embryonic samples.
#' 
#' Two time series are included in this dataset.
#' 
#'  - `H`: \insertCite{hashimshony2015spatiotemporal}{wormRef}
#'  - `L`: \insertCite{levin2016mid}{wormRef}
#' 
#' 
#' The gene expression matrix is \eqn{log(X + 1)} of quantile-normalized raw data.
#' 
#' @docType data
#' 
#' @eval data_format()
#' 
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50548}{GSE50548} 
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60755}{GSE60755}
#' 
#' @references
#'    \insertAllCited{}
#' 
#' 
#' @importFrom Rdpack reprompt
"Cel_embryo"
