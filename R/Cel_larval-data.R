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
#' @section Format :
#'  A list with :
#'  - `$g`: The gene expression matrix (genes as rows, samples as columns)
#'  - `$p`: A dataframe of phenotypic data on the samples, usually :
#'     - `sname` : The sample names (and colnames of the gene expression matrix).
#'     - `age` : Sample developmental age. May be scaled to join multiple datasets.
#'     - `cov` : A factor separating the samples into groups if needed (*e.g* multiple datasets, batches); else NULL.
#'     - `age_ini` : Sample chronological age, as given in the literature.
#'     - `accession` : Sample accession ID for GEO/ArrayExpress.
#'   - `$df`: The degree of freedom to use for the PLSR interpolation model.
#' 
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