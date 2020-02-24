sapply(c("curl", "limma", "Biobase", "utils", "RAPToR", "stats", "ica"), 
       requireNamespace, quietly = T) 

# utils for id/format conversion
load("data/Cel_genes.RData")

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

ctrl_probes <- which(grepl('Corner', RG$genes$SystematicName))
Xa <- stats::aggregate(X[-ctrl_probes, ], by = list(RG$genes$SystematicName[-ctrl_probes]), mean)
X <- Xa[, -1]
rownames(X) <- as.character(Xa[,1])

X <- RAPToR::format_ids(X, Cel_genes, from = "sequence_name", to = "wb_id")
rm(Xa, RG)

# pheno data
P <- P[P$Characteristics.age. != "not available", ]
P <- P[, c("Source.Name", "Characteristics.strain.", 
           "Characteristics.age.", "Characteristics.infect.")]
colnames(P)  <- c("sname", "strain", "age_ini", "infect")
P$infect <- factor(P$infect, levels = c("Orsay virus", "none"), labels = c("I", "NI"))
P$cov <- factor(paste0(P$strain, '.', P$infect))
P$age_ini <- as.numeric(P$age_ini)

P <- P[order(P$age_ini),]
X <- X[,P$sname]

### cleanup
if(dir.exists(raw_dir)){
  file.remove(unlist(raw_files))
  file.remove(rawzips)
  file.remove(raw_dir)
}

rm(raw_files, rawzips, 
   p_url, ctrl_probes, 
   to_dl, to_read, 
   Cel_genes, raw_dir)
  

### build Cel_YA_1

X <- limma::normalizeBetweenArrays(X, method = "quantile")
X <- log(X + 1)

# stage early N2 samples
load("data/Cel_larval.RData")
m_larv <- RAPToR::ge_im(X = Cel_larval$g, p = Cel_larval$p, formula = Cel_larval$geim_params$formula,
                        method = Cel_larval$geim_params$method, dim_red = Cel_larval$geim_params$dim_red,
                        nc = Cel_larval$geim_params$nc)

ndat <- data.frame(age = seq(min(Cel_larval$p$age), max(Cel_larval$p$age), l = 500),
                   cov = rep(Cel_larval$p$cov[1], 500))

r_larv <- list(interpGE = predict(m_larv, ndat), time.series = ndat$age)

sN2 <- P$strain == "N2"
to_stage <- sN2 & P$age_ini < max(r_larv$time.series)

ae_young_N2 <- RAPToR::ae(X[,to_stage], r_larv$interpGE, r_larv$time.series,
                          nb.cores = 3)

# adjust the age of the N2 samples
dat <- cbind(P[to_stage,], aeN2 = ae_young_N2$age.estimates[,1])
lm_N2 <- lm(aeN2 ~ age_ini, data = dat)

P$age <- predict(lm_N2, P)
P$age[to_stage] <- ae_young_N2$age.estimates[,1]

# build temp N2 reference and stage all samples
pca <- stats::prcomp(X[, sN2], rank = 20)
nc <- sum(summary(pca)$importance[3, ] < .999) + 1

m <- ge_im(X[, sel], P[sel,], formula = "X ~ s(age, bs = 'cr')", dim_red = "ica", nc = nc)
ndat <- data.frame(age = seq(min(P[sel,"age"]), max(P[sel,"age"]), l = 500))

rN2 <- list(g = predict(m, ndat), ts = ndat$age)
ae_N2 <-  RAPToR::ae(X, rN2$g, rN2$ts, nb.cores = 3)

# par(mfrow = c(2,2))
# plot(P$age_ini, ae_N2$age.estimates[,1],
#      main = "Chron. vs ae", xlab = "Age", ylab = "ae",
#      col = P$strain, lwd = 2)
# legend("bottomright", bty = 'n', lwd = 3, col = 1:3, legend = levels(P$strain),
#        lty = NA, pch = 1, text.font = 2)
# 
# invisible(sapply(levels(P$strain), function(l){
#    s <- which(P$strain == l)
#    plot(P$age_ini[s], ae_N2$age.estimates[s,1],
#         main = paste0("Chron. vs ae (", l, ")"), xlab = "Age", ylab = "ae",
#         col = which(l == levels(P$strain)), lwd = 2)
# }))

P$age[!sN2] <- ae_N2$age.estimates[!sN2, 1]

P$accession <- "E-MTAB-7573"
P <- P[, c("sname", "age", "cov", "age_ini", "accession")]
X <- X[, P$sname]


# Get nc for final reference building
pca <- stats::prcomp(X, rank = 45)
nc <- sum(summary(pca)$importance[3, ] < .999) + 1


Cel_YA_1 <- list(g = X,
                 p = P,
                 geim_params = list(formula = "X ~ s(age, bs = 'cr') + cov",
                                    method = "gam",
                                    dim_red = "ica",
                                    nc = nc)
)

# save object to data
save('Cel_YA_1', file = "data/Cel_YA_1.RData", compress = "xz")
rm(X, P,
   sN2, to_stage,
   ae_young_N2, ae_N2, 
   rN2, lm_N2, 
   r_larv, dat,
   Cel_larval, Cel_YA_1)
