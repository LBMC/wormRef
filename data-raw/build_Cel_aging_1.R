sapply(c("GEOquery", "limma", "Biobase", "utils", "RAPToR", "stats"), 
       requireNamespace, quietly = T) 

# utils for id/format conversion
load("data/Cel_genes.RData")
load("data/Cel_larval.RData")
source("data-raw/convert2tpm.R")


data_folder <- "data-raw/"

geo_id <- "GSE93826"
geo_obj <- GEOquery::getGEO(geo_id)[[1]]

# get expression data
sfs <- GEOquery::getGEOSuppFiles(geo_id, fetch_files = F, makeDirectory = F)

tmpfolder <- file.path(data_folder, "byrne2020tmp")
dir.create(tmpfolder)
destfile <- file.path(tmpfolder, sfs$fname[1])
download.file(sfs$url[1], destfile = destfile)

untar(destfile, exdir = tmpfolder)

flist <- list.files(tmpfolder)
g <- lapply(seq_along(flist)[-1], function(i){
  gi <- read.table(gzfile(file.path(tmpfolder,flist[i])), 
                   h=F, row.names = 1)
  colnames(gi) <- gsub("(GSM246\\d+)_.*", "\\1", flist[i])
  return(gi)
})

g <- do.call(cbind, g)

# format and convert to tpm
g <- RAPToR::format_ids(g, wormRef::Cel_genes, from = "wb_id", 
                        to = "wb_id", aggr.fun = sum)

g <- raw2tpm(g, wormRef::Cel_genes$transcript_length[
  match(rownames(g), wormRef::Cel_genes$wb_id)])
g <- g[apply(g, 1, sum) > 0, ] # keep expressed genes only
g <- g[-3123, ] # remove 0-variance artefact gene

# get pheno data
p <-  Biobase::pData(geo_obj)
p <- p[, c("title", "geo_accession", "age:ch1")]
colnames(p) <- c("sname", "accession", "age")
p$age <- as.numeric(as.character(
  factor(p$age, levels = paste0('day ', c(3,6,9,12,15,18)),
         labels = c(3,6,9,12,15,18))
))

p <- p[order(p$age), ]
g <- g[, p$accession]



# normalize and log
g <- log1p(limma::normalizeBetweenArrays(g, method = "quantile"))

# compute correlation to keep strong aging signal genes
gn_cor <- apply(g, 1, cor, y=p$age, method = "spearman")
selg <- (abs(gn_cor)>sqrt(1/3))


Cel_aging_1 <- list(g = g[selg,],
                    p = p,
                    full.g = g, # keep also full expr. data
                    geim_params = list(formula = "X ~ s(age, bs = 'cr', k=4)",
                                       method = "gam",
                                       dim_red = "pca",
                                       nc = 1),
                    t.unit = "days of adulthood (20C)",
                    cov.levels = NULL,
                    metadata = list("organism" = "C. elegans",
                                    "profiling" = "whole-organism, bulk RNAseq",
                                    "source" = "GSE93826",
                                    "condition" = "rrf-3(pk1436), liquid culture")
)

# save object to data
save('Cel_aging_1', file = "data/Cel_aging_1.RData", compress = "xz")


# cleanup
file.remove(file.path(tmpfolder, flist))
file.remove(tmpfolder, destfile)
rm(g, p, tmpfolder, flist, geo_id, geo_obj, sfs, destfile,
   gn_cor, selg,
   Cel_genes, Cel_aging_1)

