# wormRef

This repo is an R package with the nematode development references for age estimation using the [`RAPToR` tool](https://github.com/LBMC/RAPToR).

To install the package, run the following in your R console :

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("LBMC/wormRef", build_vignettes = T)
```

<hr>

## Update info
### v0.5
 - Built a function factory to simplify ref. building code
 - Added reference metadata to comply with RAPToR v1.2 update.
 - Added time range/unit info in list_refs 
 - Added an aging reference (`Cel_aging_1`)
 
### v0.4
 - Added Meeuse et al. (2020) reference spanning L1 to adult development

### v0.3
#### v0.3.1
 - Fixed plot_refs not diplaying
 - Removed doc typos
 
#### v0.3.0
 - Updated reference parameters with RAPToR bug fix
 - Changed all RPKM to TPM
 - Updated the embryo reference with stricter quality filter on samples 
 
### v0.2
 - Updated the internal structure to match the new `RAPToR` interface witth data-packages.
 - Re-evaluated optimal interpolation of reference series using the new GEIM approach of `RAPToR`.
 
### v0.1
 - Created the data package. 
