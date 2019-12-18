
lapply(c("Cel_embryo", "Cel_larval", "Cel_YA_1", "Cel_YA_2", "Cel_genes"), function(d){
  d_path <- paste0("data/",d,".RData")
  load(d_path)
  save(list = d, file = d_path, compress = "xz")
  rm(list = d)
})

gc()
