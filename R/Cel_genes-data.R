#' Cel_genes
#' 
#' A gene ID table for *C. elegans*.
#' 
#' This table was built directly from the paraSite biomaRt.
#' 
#' 
#' @docType data
#' 
#' @section Format:
#' A dataframe with the folowing columns :
#' 
#'  - `wb_id`: WBGene ID (*e.g* `WBGene00019636`)
#'  - `transcript_name`: Transcript name (*e.g* `K10F12.4a.1`)
#'  - `sequence_name`: Sequence name (*e.g* `K10F12.4`) - this is most frequently used in transcriptomic datasets
#'  - `public_name`: Gene name (*e.g* `gsto-3`)
#'  - `transcript_length`: Transcript length from 5'UTR start to 3'UTR end (*e.g* `3161`)
#' 
#' 
#' @source \href{https://parasite.wormbase.org/info/Tools/biomart.html}{paraSite biomaRt} 
#' 
"Cel_genes"
