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
  - `$df`: The degree of freedom to use for the PLSR interpolation model.
  - `$nc`: The number of PLS components to use for the prediction of the PLSR interpolation model
"
  return(strsplit(txt, split = "\n")[[1]])
}