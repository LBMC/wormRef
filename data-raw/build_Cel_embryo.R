sapply(c("GEOquery", "limma", "Biobase", "utils", "RAPToR"), 
       requireNamespace, quietly = T) 

# utils for id/format conversion
load("data/Cel_genes.RData")
# source("data-raw/raw2rpkm.R")
source("data-raw/convert2tpm.R")


### get Levin data
geo_id_L <- "GSE60755"

g_url_L <- GEOquery::getGEOSuppFiles(geo_id_L, makeDirectory = FALSE, fetch_files = FALSE)
g_file_L <- "data-raw/levin.txt.gz"
utils::download.file(url = as.character(g_url_L$url[1]), destfile = g_file_L)

X_L <- read.table(gzfile(g_file_L), header = T, row.names = 1, sep = "\t")

# convert to tpm & wb_id
X_L <- X_L[rownames(X_L)%in%Cel_genes$sequence_name,]
X_L <- raw2tpm(rawcounts = X_L, genelengths = Cel_genes$transcript_length[match(rownames(X_L), Cel_genes$sequence_name)])

X_L <- RAPToR::format_ids(X_L, Cel_genes, from = "sequence_name", to = "wb_id")


P_L <- GEOquery::getGEO(geo_id_L)[[1]]

P_L <- Biobase::pData(P_L)
P_L <- P_L[, c("title", "geo_accession", "time point (minutes after 4-cell):ch1")]
colnames(P_L) <- c("sname", "accession", "age")
P_L$age <- as.numeric(P_L$age)
P_L$age_ini <- P_L$age


# formatting
colnames(X_L) <- gsub("Metazome_CE_timecourse_", "", colnames(X_L))
P_L$sname <- gsub("Metazome_CE_timecourse_", "", P_L$sname)

P_L <- P_L[order(P_L$age), ]
X_L <- X_L[, P_L$sname]


### cleanup
if(file.exists(g_file_L))
   file.remove(g_file_L)
rm(fpkm2tpm, raw2tpm, g_url_L, g_file_L, geo_id_L)



### build Cel_embryo object
# filter bad Levin samples
f_lev <- c(which(0.67 > apply(RAPToR::cor.gene_expr(X_L, X_L), 1, quantile, probs = .99)),
           which(P_L$sname %in% c("sample_0097", "sample_0099", "sample_0046", "sample_0053", 
                                  "sample_0054", "sample_0073", "sample_0100", "sample_0055", 
                                  "sample_0056", "sample_0010", "sample_0105", "sample_0009", 
                                  "sample_0106", "sample_0058", "sample_0059", "sample_0001", 
                                  "sample_0002", "sample_0003", "sample_0004", "sample_0125")))
# Named samples are removed bc they are clear outliers from PCA components, w.r.t dynamics.


P_L <- P_L[-f_lev, c("sname", "age", "accession", "age_ini")]
X_L <- X_L[, -f_lev]


# Normalize
X_L <- limma::normalizeBetweenArrays(X_L, method = "quantile")
X_L <- log1p(X_L)





# get nc for final reference building
tXc <- scale(t(X_L), center = TRUE, scale = FALSE)
pca <- summary(stats::prcomp(tXc, rank = 25))
nc <- sum(pca$importance[3, ] < .7) + 1 # .90 bc of medium quality

# ks <- c(seq(4,16, 2))
# flist <- as.list(c(paste0("X ~ s(age, bs = 'cr', k=",ks,")"), "X ~ s(age, bs = 'cr')"))
# cv_L <- ge_imCV(X = X_L, p = P_L, formula_list = flist, method = "gam", dim_red = "pca", nc = nc)
# 
# plot(cv_L, names = paste0("k=", c(ks, 'n')))

Cel_embryo <- list(g = X_L,
                   p = P_L,
                   geim_params = list(formula = "X ~ s(age, bs = 'cr', k = 9)",
                                      method = "gam",
                                      dim_red = "pca",
                                      nc = nc),
                   t.unit = "min past 4-cell (20C)",
                   cov.levels = NULL,
                   metadata = list("organism" = "C. elegans",
                                   "profiling" = "single-embryo RNAseq",
                                   "source" = "GSE60755")
)


# save object to data
save('Cel_embryo', file = "data/Cel_embryo.RData", compress = "xz")
rm(X_L, P_L, tXc,
   f_lev, pca, nc,
   Cel_genes, Cel_embryo)
