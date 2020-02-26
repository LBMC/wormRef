# wormRef

This repo is an R package with the nematode development references for age estimation using the [`RAPToR` tool](https://github.com/LBMC/RAPToR).

To install the package, you can use the `devtools` R package. This should be done in your R console :

```r
library(devtools)
devtools::install_github("LBMC/wormRef")
```

If you don't have `devtools` installed, you can do the following :
```r
install.packages("devtools")
```

<hr>

## Update info

### v0.2
 - Updated the internal structure to match the new `RAPToR` interface witth data-packages.
 - Re-evaluated optimal interpolation of reference series using the new GEIM approach of `RAPToR`.
 
### v0.1.4
 - Created the data package. 
