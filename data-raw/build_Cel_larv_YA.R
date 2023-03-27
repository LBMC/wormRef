sapply(c("GEOquery", "limma", "Biobase", "utils", "RAPToR", "stats"), 
       requireNamespace, quietly = T) 

# utils for id/format conversion
load("data/Cel_genes.RData")
load("data/Cel_larval.RData")
source("data-raw/convert2tpm.R")

# Get Meeuse2020 data
geo_id <- "GSE130811"
g_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130811&format=file&file=GSE130811%5Fexpr%5FmRNA%5FCE10%5Fcoding%2Etab%2Egz"
g_f <- "data-raw/mee2020.tab.gz"

utils::download.file(g_url, destfile = g_f)

g <- read.table(gzfile(g_f), h=T, row.names = 1, as.is = T, sep = "\t")


# format ids
g <- RAPToR::format_ids(g, Cel_genes, from = "wb_id", to = "wb_id")

# convert to tpm
g <- raw2tpm(g[,-1], genelengths = Cel_genes[match(rownames(g), Cel_genes$wb_id), "transcript_length"])
# remove strange duplicate samples
g <- g[,1:44]


# pheno Data
p <- Biobase::pData(GEOquery::getGEO(geo_id)[[1]])
p <- p[, c("title", "geo_accession", "developmental stage:ch1")]
colnames(p) <- c("sname", "accession", "age_ini")


p$age_ini <- as.numeric(gsub("(\\d+)h", "\\1", p$age_ini) )
p$sname <- make.names(gsub("RNAseq N2 (\\d+h)", "\\1r", p$sname))

#all(p$sname==colnames(g))



# normalize & log
g <- log1p(limma::normalizeBetweenArrays(g, method = "quantile"))


# Estimate developmental speed diff. with 20C larval reference to adjust age and match 20C development
m_larv <- RAPToR::ge_im(X = Cel_larval$g, p = Cel_larval$p, formula = Cel_larval$geim_params$formula,
                        method = Cel_larval$geim_params$method, dim_red = Cel_larval$geim_params$dim_red,
                        nc = Cel_larval$geim_params$nc)
ndat <- data.frame(age = seq(min(Cel_larval$p$age), max(Cel_larval$p$age), l = 500),
                   cov = rep(Cel_larval$p$cov[60], 500))

r_larv <- list(interpGE = predict(m_larv, ndat), time.series = ndat$age)

sel <- p$age_ini*1.5 < 48
ae_young <- RAPToR::ae(g[, sel], r_larv$interpGE, r_larv$time.series, prior = p$age_ini[sel]*1.5, prior.params = 10)



lm_y <- stats::lm(ae_young$age.estimates[,1]~ 0 + p$age_ini[sel])
p$age <- lm_y$coefficients[1] * p$age_ini




# get nc for ge_im
pca <- summary(stats::prcomp(t(g), scale= F, center = T, rank = 40))
nc <- sum(pca$importance["Cumulative Proportion",] < .99) + 1


p <- p[, c("sname", "age", "age_ini", "accession")]

Cel_larv_YA <- list(g = g,
                   p = p,
                   geim_params = list(formula = "X ~ s(age, k = 25, bs = 'cr')",
                                      method = "gam",
                                      dim_red = "ica",
                                      nc = nc),
                   t.unit = "h past egg-laying (20C)",
                   cov.levels = NULL,
                   metadata = list("organism" = "C. elegans",
                                   "profiling" = "whole-organism, bulk RNAseq",
                                   "source" = "GSE130811")
)

# save object to data
save('Cel_larv_YA', file = "data/Cel_larv_YA.RData", compress = "xz")


file.remove(g_f)
rm(g, p, g_f, g_url, sel, geo_id,
   m_larv, ndat, r_larv, ae_young,
   nc, pca, lm_y,
   Cel_genes, Cel_larval, Cel_larv_YA)
