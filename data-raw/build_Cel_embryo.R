sapply(c("GEOquery", "limma", "Biobase", "utils", "RAPToR"), 
       requireNamespace, quietly = T) 

# utils for id/format conversion
load("data/Cel_genes.RData")
source("data-raw/raw2rpkm.R")



### get Hashimshony data
geo_id_H <- "GSE50548"

g_url_H <- GEOquery::getGEOSuppFiles(geo_id_H, makeDirectory = FALSE, fetch_files = FALSE)
g_file_H <- "data-raw/hash.tab.gz"
utils::download.file(url = as.character(g_url_H$url[2]), destfile = g_file_H)

X_H <- read.table(gzfile(g_file_H), h=T, sep = '\t', stringsAsFactors = F, row.names = 1)

# convert to rpkm & wb_id
X_H <- RAPToR::format_ids(X_H, Cel_genes, from = "sequence_name", to = "wb_id")
X_H <- raw2rpkm(X = X_H, gene.length = Cel_genes, id.col = "wb_id", l.col = "transcript_length")


# pheno data
P_H <- GEOquery::getGEO('GSE50548', getGPL = F, GSEMatrix = F)
P_H <- GEOquery::GSMList(P_H)
P_H <- lapply(P_H, GEOquery::Meta)

P_H <- do.call('rbind', lapply(P_H, function(p){
  unlist(p[c("title", "source_name_ch1", "geo_accession", "description", "characteristics_ch1")])
}))

# keep only full embryo samples
P_H <- P_H[grepl('embryo', P_H[,"source_name_ch1"]),]
P_H <- P_H[, c("title", "geo_accession", "description")]

P_H <- as.data.frame(P_H, stringsAsFactors = F)

P_H$age <- as.numeric(gsub('(\\d+)\\s*-*.*', '\\1', as.character(P_H$description)))
P_H$age[1:2] <- c(-50,-30)

P_H$title[1:2] <- paste0('X', P_H$title[1:2])
P_H$description <- NULL

X_H <- X_H[, P_H$title]



### get Levin data
geo_id_L <- "GSE60755"

g_url_L <- GEOquery::getGEOSuppFiles(geo_id_L, makeDirectory = FALSE, fetch_files = FALSE)
g_file_L <- "data-raw/levin.txt.gz"
utils::download.file(url = as.character(g_url_L$url[1]), destfile = g_file_L)

X_L <- read.table(gzfile(g_file_L), header = T, row.names = 1, sep = "\t")

# convert to rpkm & wb_id
X_L <- RAPToR::format_ids(X_L, Cel_genes, from = "sequence_name", to = "wb_id")
X_L <- raw2rpkm(X = X_L, gene.length = Cel_genes, id.col = "wb_id", l.col = "transcript_length")


# pheno data
P_L <- GEOquery::getGEO(geo_id_L)[[1]]

P_L <- Biobase::pData(P_L)
P_L <- P_L[, c("title", "geo_accession", "time point (minutes after 4-cell):ch1")]
colnames(P_L)[3] <- "age"
P_L$age <- as.numeric(P_L$age)

# formatting
colnames(X_L) <- gsub("Metazome_CE_timecourse_", "", colnames(X_L))
P_L$title <- gsub("Metazome_CE_timecourse_", "", P_L$title)

P_L <- P_L[order(P_L$age), ]
X_L <- X_L[, P_L$title]


### cleanup
if(file.exists(g_file_H))
  file.remove(g_file_H)
if(file.exists(g_file_L))
  file.remove(g_file_L)
rm(raw2rpkm, g_url_H, g_url_L, g_file_H, g_file_L, geo_id_H, geo_id_L)



### build Cel_embryo object

# filter bad Levin samples
ccl <- RAPToR::cor.gene_expr(X_L, X_L)
f_lev <- which(0.6 > apply(ccl, 1, quantile, probs = .99))

# join datasets
X <- RAPToR::format_to_ref(X_H, X_L[, -f_lev])
X <- cbind(X[[1]], X[[2]])

X <- limma::normalizeBetweenArrays(X, method = "quantile")
X <- log(X + 1)


P_H$cov <- "H"
P_L$cov <- "L"

P <- rbind(P_H, P_L[-f_lev, ])

# formatting
P$cov <- factor(P$cov, levels = c("H", "L"))
P$age_ini <- P$age
colnames(P) <- c("sname", "accession", "age", "cov", "age_ini")
P <- P[, c("sname", "age", "cov", "age_ini", "accession")]
X <- X[, P$sname]

Cel_embryo <- list(g = X,
                   p = P,
                   df = 19)

# save object to data
save('Cel_embryo', file = "data/Cel_embryo.RData")
rm(X_H, X_L, X, 
   P_H, P_L, P,
   ccl, f_lev,
   Cel_genes, Cel_embryo)
