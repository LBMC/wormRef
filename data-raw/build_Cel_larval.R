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





### get Hendriks data
geo_id_H <- "GSE52861"

g_url_H <- GEOquery::getGEOSuppFiles(geo_id_H, makeDirectory = FALSE, fetch_files = FALSE)
g_file_H <- "data-raw/hendriks.txt.gz"
utils::download.file(url = as.character(g_url_H$url[2]), destfile = g_file_H)

X_H <- read.table(gzfile(g_file_H), h=T, sep = '\t', stringsAsFactors = F, row.names = 1)

# convert to rpkm & wb_id
X_H <- wormAge::format_ids(X_H, Cel_genes, from = "wb_id", to = "wb_id")
X_H <- raw2rpkm(X = X_H, gene.length = Cel_genes, id.col = "wb_id", l.col = "transcript_length")


# pheno data
P_H <- Biobase::pData(GEOquery::getGEO(geo_id_H, getGPL = F)[[1]])

# filter relevant fields/samples
P_H <- P_H[(P_H$`strain:ch1` == 'N2') & (P_H$`growth protocol:ch1` == 'Continuous'), ]
P_H <- P_H[, c("title", "geo_accession", "time in development:ch1")]

# get age 
P_H$age <- as.numeric(sub('(\\d+)\\shours', '\\1', P_H$`time in development:ch1`))

# formatting
P_H$title <- gsub('RNASeq_polyA_', '', 
                  gsub('hr', 'h', 
                       gsub('-', '.', fixed = T, as.character(P_H$title))))
colnames(X_H) <- gsub('RNASeq_polyA_','', colnames(X_H))

P_H$cov <- "H"
X_H <- X_H[, P_H$title]



### cleanup
if(file.exists(g_file_O))
  file.remove(g_file_O)
if(file.exists(g_file_H))
  file.remove(g_file_H)
rm(raw2rpkm, g_url_O, g_url_H, g_file_O, g_file_H, geo_id_O, geo_id_H)



### build Cel_larval

# join datasets
X <- wormAge::format_to_ref(X_O, X_H)
X <- cbind(X[[1]], X[[2]])

X <- limma::normalizeBetweenArrays(X, method = "quantile")
X <- log(X + 1)

P <- rbind(P_O[, c("title", "geo_accession", "age", "cov")], P_H[, c("title", "geo_accession", "age", "cov")])

# formatting
P$cov <- factor(P$cov, levels = c("O.20", "O.25", "H"))
P$age_ini <- P$age
