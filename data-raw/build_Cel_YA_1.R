sapply(c("curl", "limma", "Biobase", "utils", "wormAge", "stats"), 
       requireNamespace, quietly = T) 


### Get Sterken data
p_url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7573/E-MTAB-7573.sdrf.txt"
P <- read.table(p_url, sep = '\t', h = T, quote = "\"", comment.char = "", stringsAsFactors = F)

raw_dir <- "data-raw/sterken_raw/"
dir.create(raw_dir, showWarnings = FALSE)
to_dl <- unique(P$Comment..ArrayExpress.FTP.file.)

rawzips <- sapply(to_dl, function(f){
  curl::curl_download(f, destfile = paste0(raw_dir, strsplit(f, 'E-MTAB-\\d+/')[[1]][2]))
})

raw_files <- sapply(rawzips, utils::unzip, exdir = raw_dir)
to_read <- unique(P$Array.Data.File)

RG <- limma::read.maimages(files = to_read, path = raw_dir, source = 'agilent')
G <- RG$G
R <- RG$R

Gi <- match(colnames(G), P[P$Label == 'Cy3', "Assay.Name"])
Ri <- match(colnames(R), P[P$Label == 'Cy5', "Assay.Name"])

colnames(G) <- (P[P$Label == 'Cy3',])[Gi, "Source.Name"]
colnames(R) <- (P[P$Label == 'Cy5',])[Ri, "Source.Name"]

X <- cbind(G, R)

rm(R, G, Ri, Gi)

