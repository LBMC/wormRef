sapply(c("GEOquery", "limma", "Biobase", "utils", "RAPToR", "stats"), 
       requireNamespace, quietly = T) 

# utils for id/format conversion
load("data/Cel_genes.RData")

### get Reinke data
geo_ids <- c("GSE726", "GSE727", "GSE728", "GSE729", 
             "GSE730", "GSE731", "GSE732", "GSE733",
             "GSE734", "GSE735", "GSE736", "GSE737")

Gs <- lapply(geo_ids, GEOquery::getGEO, GSEMatrix = F)

# pheno data
P <- do.call(rbind, lapply(Gs, function(g){
  data.frame(
      apply(do.call(rbind, lapply(GEOquery::GSMList(g), GEOquery::Meta))[, c("geo_accession", "source_name_ch2")],
            2, unlist),
    stringsAsFactors = F)
}))

# time is +3 hours per timepoint (cf. paper)
P$age_ini <- (as.numeric(gsub('TP', '', gsub('(TP\\d+).*', '\\1', P$source_name_ch2))) - 1) * 3
P$cov <- factor(gsub('.*prep\\s(\\d).*', '\\1', P$source_name_ch2))



# geno data
X <- do.call(cbind, lapply(Gs, function(g){
  do.call(cbind, lapply(GEOquery::GSMList(g), function(gl){
    GEOquery::Table(gl)[, "CH2_median"] #- GEOquery::Table(gl)[, "CH2_back"]
  }))
}))

X[X<0] <- 0


# get/convert probe ids
probe_ids <- as.character(GEOquery::Table(GEOquery::GSMList(Gs[[1]])[[1]])[,"ID_REF"])
gpl <- GEOquery::Table(GEOquery::GPLList(Gs[[1]])[[1]])
rownames(gpl) <- as.character(gpl$ID)

rownames(X) <- gpl[probe_ids, "ORF"]

X <- RAPToR::format_ids(X, Cel_genes, from = "sequence_name", to = "wb_id")

X <- X[, P$geo_accession]

# log/normalize
X <- limma::normalizeBetweenArrays(X, method = "quantile")
X <- log1p(X)



# adjust ages to 20C development from hatching
load("data/Cel_larval.RData")
m_larv <- RAPToR::ge_im(X = Cel_larval$g, p = Cel_larval$p, formula = Cel_larval$geim_params$formula,
                        method = Cel_larval$geim_params$method, dim_red = Cel_larval$geim_params$dim_red,
                        nc = Cel_larval$geim_params$nc)

ndat <- data.frame(age = seq(min(Cel_larval$p$age), max(Cel_larval$p$age), l = 500),
                   cov = rep(Cel_larval$p$cov[1], 500))

r_larv <- list(interpGE = predict(m_larv, ndat), time.series = ndat$age)

to_stage <- (35 + P$age_ini * 1.5) < max(r_larv$time.series)
ae_young <- RAPToR::ae(X[,to_stage], r_larv$interpGE, r_larv$time.series,
                       bootstrap.n = 1)

dat <- cbind(P[to_stage,], ae = ae_young$age.estimates[,1])
lm_r <- lm(ae ~ age_ini + cov, data = dat)

P$age <- predict(lm_r, P)
P$age[to_stage] <- ae_young$age.estimates[,1]


# formatting
P$accession <- P$geo_accession
colnames(P) <- c("sname", "n", "age_ini", "cov", "age", "accession")

P <- P[, c("sname", "age", "cov", "age_ini", "accession")]
X <- X[, P$sname]


# get nc for final reference building
tXc <- scale(t(X), scale = FALSE, center = TRUE)
pca <- summary(stats::prcomp(tXc, rank = 20))
nc <- sum(pca$importance[3, ] < .8) + 1 # only .9 bc of bad quality


Cel_YA_2 <- list(g = X,
                 p = P,
                 geim_params = list(formula = "X ~ s(age, bs = 'tp') + cov",
                                    method = "gam",
                                    dim_red = "pca",
                                    nc = nc)
)

# save object to data
save('Cel_YA_2', file = "data/Cel_YA_2.RData", compress = "xz")

rm(X, P, tXc, 
   Gs, gpl, geo_ids,
   probe_ids, to_stage,
   ae_young, dat, lm_r, 
   Cel_larval, r_larv, m_larv, 
   ndat, pca, nc,
   Cel_YA_2, Cel_genes)

