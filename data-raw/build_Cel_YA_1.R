sapply(c("curl", "limma", "Biobase", "utils", "RAPToR", "stats"), 
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
r_larv <- RAPToR::plsr_interpol(Cel_larval$g, Cel_larval$p$age, 
                                df = Cel_larval$df, covar = Cel_larval$p$cov, 
                                topred = "O.20", n.inter = 500)

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
rN2 <- RAPToR::plsr_interpol(X[, sN2], P$age_ini[sN2], df = 5, covar = P$infect[sN2], 
                             topred = 'NI', n.inter = 200)
ae_N2 <-  RAPToR::ae(X, rN2$interpGE, rN2$time.series, nb.cores = 3)

P$age[!sN2] <- ae_N2$age.estimates[!sN2, 1]

P$accession <- "E-MTAB-7573"
P <- P[, c("sname", "age", "cov", "age_ini", "accession")]
X <- X[, P$sname]

Cel_YA_1 <- list(g = X,
                 p = P,
                 df = 5,
                 nc = 8)

# save object to data
save('Cel_YA_1', file = "data/Cel_YA_1.RData")
rm(X, P,
   sN2, to_stage,
   ae_young_N2, ae_N2, 
   rN2, lm_N2, 
   r_larv, dat,
   Cel_larval, Cel_YA_1)
