sapply(c("GEOquery", "limma", "Biobase", "utils", "wormAge"), 
       requireNamespace, quietly = T) 

# utils for id/format conversion
load("data/Cel_genes.RData")
source("data-raw/raw2rpkm.R")

### get Oudenaarden data
geo_id_O <- "GSE49043"

g_url_O <- GEOquery::getGEOSuppFiles(geo_id_O, makeDirectory = FALSE, fetch_files = FALSE)
g_file_O <- "data-raw/ouden.txt.gz"
utils::download.file(url = as.character(g_url_O$url[1]), destfile = g_file_O)

X_O <- read.table(gzfile(g_file_O), h=T, sep = '\t', stringsAsFactors = F)
X_O <- X_O[!duplicated(X_O$GENEID), ] # remove duplicate rows
rownames(X_O) <- X_O$GENEID
X_O <- X_O[, -c(1:2)]

# format ids
X_O <- wormAge::format_ids(X_O, Cel_genes, from = "sequence_name", to = "wb_id")


# pheno Data
P_O <- Biobase::pData(GEOquery::getGEO(geo_id_O, getGPL = F)[[1]])

# filter relevant fields/samples
P_O <- P_O[, c("title", "geo_accession", "strain:ch1", "temperature:ch1")]
P_O <- P_O[P_O$`strain:ch1`=='N2',]
P_O$`strain:ch1` <- NULL
P_O$title <-  as.character(P_O$title)

# get age from sample names
P_O$age <- as.numeric(sub('DH\\d_N2_(\\d+)_?(\\d)?', '\\1\\.\\2', P_O$title))

P_O$cov <- paste0("O.", gsub('C', '', P_O$`temperature:ch1`))

X_O <- X_O[, P_O$title]

