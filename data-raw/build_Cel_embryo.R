sapply(c("GEOquery", "limma", "Biobase", "utils", "wormAge"), 
       requireNamespace, quietly = T) 

# Get Hashimshony data
geo_id_H <- "GSE50548"

g_url <- GEOquery::getGEOSuppFiles(geo_id_H, makeDirectory = FALSE, fetch_files = FALSE)
g_file <- "data-raw/hash.tab.gz"
utils::download.file(url = as.character(g_url$url[2]), destfile = g_file)

X_H <- read.table(gzfile(g_file), h=T, sep = '\t', stringsAsFactors = F, row.names = 1)

# convert to rpkm & wb_id
load("data/Cel_genes.RData")
source("data-raw/raw2rpkm.R")

X_H <- X_H[rownames(X_H)%in%Cel_genes$transcript_name, ]
X_H <- raw2rpkm(X = X_H, gene.length = Cel_genes, id.col = "transcript_name", l.col = "transcript_length")
X_H <- wormAge::format_ids(X_H, Cel_genes, from = "transcript_name", to = "wb_id")