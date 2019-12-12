datasets_s <- "name, tstart, tend, tunit
Cel_embryo, -50, 840, min past 4C
Cel_larval, 0, 55, hours post-hatching
Cel_YA_1, 39, 72, hours post-hatching
Cel_YA_2, 48, 83, hours post-hatching
"
devstages_s <-"devstage, tstart, tend, tunit
1C, -50, -20, min past 4C
2C, -20,   0, min past 4C
4C,   0,  20, min past 4C
Gastrulation, 150, 350, min past 4C
Comma, 400, 460, min past 4C
2-fold, 460, 520, min past 4C
4-fold, 520, 600, min past 4C
L1, 0,  20, hours post-hatching
L2, 20, 30, hours post-hatching
L3, 30, 38, hours post-hatching
L4, 38, 50, hours post-hatching
Young Adult, 50, 65, hours post-hatching
Adult, 65, 80, hours post-hatching
"
devst <- read.table(text = devstages_s, h = T, sep = ',', stringsAsFactors = F)
devst$tunit <- factor(devst$tunit)

dss <- read.table(text = datasets_s, h = T, sep = ',', stringsAsFactors = F)
dss$tunit <- factor(dss$tunit)

Cel_devstages <- list(
  devstages = devst, 
  datasets = dss)

save("Cel_devstages", file = "data/Cel_devstages.RData")

rm(devstages_s, devst,
   datasets_s, dss,
   Cel_devstages)