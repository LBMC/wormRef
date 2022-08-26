data_format <- function(){
  txt <- "
@section Format :
 A list with :
 - `$g`: The gene expression matrix (genes as rows, samples as columns)
 - `$p`: A dataframe of phenotypic data on the samples, usually :
    - `sname` : The sample names (and colnames of the gene expression matrix).
    - `age` : Sample developmental age. May be scaled to join multiple datasets.
    - `cov` : A factor separating the samples into groups if needed (*e.g* multiple datasets, batches); else NULL.
    - `age_ini` : Sample chronological age, as given in the literature.
    - `accession` : Sample accession ID for GEO/ArrayExpress.
  - `$geim_params`: A list of elements needed for the gene expression interpolation model :
    - `$formula` : the model formula.
    - `$method` : the model type (eg. gam, glm).
    - `$dim_red` : the dimension reduction method to use (pca, ica).
    - `$nc`: The number PCA/ICA components to use in the model.
  - `$t.unit`: the time series unit (e.g 'hours').
  - `$cov.levels`: a named list of potential covariate with the level to interpolate as.
  - `$metadata`: a named list of extra metadata (e.g organism, profiling technology).
"
  return(strsplit(txt, split = "\n")[[1]])
}