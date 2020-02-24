sapply(c("GEOquery", "limma", "Biobase", "utils", "RAPToR", "stats"), 
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
X_O[is.na(X_O)] <- 0

# format ids
X_O <- RAPToR::format_ids(X_O, Cel_genes, from = "sequence_name", to = "wb_id")


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
X_H <- RAPToR::format_ids(X_H, Cel_genes, from = "wb_id", to = "wb_id")
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
X <- RAPToR::format_to_ref(X_O, X_H)
X <- cbind(X[[1]], X[[2]])

X <- limma::normalizeBetweenArrays(X, method = "quantile")
X <- log(X + 1)

P <- rbind(P_O[, c("title", "geo_accession", "age", "cov")], 
           P_H[, c("title", "geo_accession", "age", "cov")])

# formatting
P$cov <- factor(P$cov, levels = c("O.20", "O.25", "H"))
P$age_ini <- P$age


# build temp 20C reference
sO20 <- P$cov == "O.20" & P$title != "DH5_N2_38" # outlier with err. time
pca <- stats::prcomp(X[, sO20], rank = 20)
summary(pca)

# select nb of components to use for interpolation
nc <- sum(summary(pca)$importance[3, ] < .999) + 1

# build geim model and predictions
m <- ge_im(X = X[, sO20], p = P[sO20,], formula = "X ~ s(age_ini, bs = 'ds')", 
           method = "gam", dim_red = "pca", nc = nc)
ndat <- data.frame(age_ini = seq(min(P[sO20, "age_ini"]), max(P[sO20, "age_ini"]), l = 500))
r20C <- list(g = predict(m, ndat), ts = ndat$age_ini)

# estimate age of samples
ae_r20 <- ae(X, r20C$g, r20C$ts, bootstrap.n = 1)

# par(mfrow = c(2,2))
# plot(P$age_ini, ae_r20$age.estimates[,1], 
#      main = "Chron. vs ae", xlab = "Age", ylab = "ae", 
#      col = P$cov, lwd = 2)
# legend("bottomright", bty = 'n', lwd = 3, col = 1:3, legend = levels(P$cov), 
#        lty = NA, pch = 1, text.font = 2)
# 
# invisible(sapply(levels(P$cov), function(l){
#    s <- which(P$cov == l)
#    plot(P$age_ini[s], ae_r20$age.estimates[s,1], 
#         main = paste0("Chron. vs ae (", l, ")"), xlab = "Age", ylab = "ae", 
#         col = which(l == levels(P$cov)), lwd = 2)
# }))

# predict age of Hendriks samples outside of rO20 range
dat <- P[P$cov == "H", ]
dat$ae <- NA
H_staged <- dat$age_ini <= (48/1.5)
dat$ae[H_staged] <- (ae_r20$age.estimates[P$cov == "H", 1])[H_staged]

lm_h <- stats::lm(ae ~ age_ini, data = dat)

dat$age <- dat$ae
dat$age[!H_staged] <- stats::predict(lm_h, dat[is.na(dat$ae),])

P$age[P$cov == "H"] <- dat$age


colnames(P) <- c("sname", "accession", "age", "cov", "age_ini")
P <- P[, c("sname", "age", "cov", "age_ini", "accession")]
X <- X[, P$sname]

# Get nc for final reference building
pca <- stats::prcomp(X, rank = 45)
nc <- sum(summary(pca)$importance[3, ] < .999) + 1


Cel_larval <- list(g = X,
                   p = P,
                   geim_params = list(formula = "X ~ s(age, bs = 'ds') + cov",
                                      method = "gam",
                                      dim_red = "pca",
                                      nc = nc)
                   )

# save object to data
save('Cel_larval', file = "data/Cel_larval.RData", compress = "xz")
rm(X_O, X_H, X, 
   P_O, P_H, P,
   dat, lm_h, 
   r20C, ae_r20,
   sO20, H_staged,
   m, ndat, pca, nc,
   Cel_genes, Cel_larval)
