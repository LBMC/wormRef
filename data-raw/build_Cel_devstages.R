

devst <- read.table(file =  file.path("data-raw", "dev_timings.csv"), h = T, sep = ',', stringsAsFactors = F)
devst$tunit <- factor(devst$tunit)

dss <- read.table(file = file.path("data-raw", "ref_timings.csv"), h = T, sep = ',', stringsAsFactors = F)
dss$tunit <- factor(dss$tunit)

Cel_devstages <- list(
  devstages = devst, 
  datasets = dss)

save("Cel_devstages", file = "data/Cel_devstages.RData")

rm(devst, dss,
   Cel_devstages)