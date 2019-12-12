#' Transcriptomic data from *C. elegans* embryonic samples
#' 
#' Two time series are included in this dataset.
#' \describe{
#'     \item{H}{Hashimshony, et al.  (2015)}
#'     \item{L}{Levin, et al.  (2016)}
#' }
#' 
#' The gene expression matrix is \eqn{log(X + 1)} of quantile-normalized raw data.
#' 
#' @docType data
#' 
#' @format A list with :
#' \describe{
#'     \item{g}{The gene expression matrix (genes as rows, samples as columns)}
#'     \item{p}{Phenotypic data on the samples, usually :
#'     \describe{
#'         \item{sname}{The sample names (and colnames of the gene expression matrix).}
#'         \item{age}{Sample developmental age. May be scaled to join multiple datasets.}
#'         \item{cov}{A factor separating the samples into groups if needed (*e.g* multiple datasets, batches); else NULL.}
#'         \item{age_ini}{Sample chronological age, as given in the literature.}
#'         \item{accession}{Sample accession ID for GEO/ArrayExpress.}
#'     }}
#'     \item{df}{The degree of freedom to use for the PLSR interpolation model.}
#' }
#' 
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50548}{GSE50548} 
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60755}{GSE60755}
#' 
#' @references 
#' \insertref{hashimshony2015spatiotemporal}{wormRef}
#' \insertref{levin2016mid}{wormRef}
"Cel_embryo"
