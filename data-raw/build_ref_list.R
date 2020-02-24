ref_list <- read.table("data-raw/ref_list.csv", h=T, sep = ',', stringsAsFactors = F)

save("ref_list", file = "data/ref_list.RData")

rm(ref_list)