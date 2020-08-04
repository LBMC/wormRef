#' Build Interpolated Gene Expression References
#' 
#' Builds the interpolation of the reference datasets.
#' **These functions are internally called by \code{\link[RAPToR]{prepare_refdata}} from RAPToR.**
#' 
#' @param n.inter the resolution of the interpolation, as in \code{seq(start, end, length.out = n.inter)}.
#' 
#' @return A list with \code{interpGE} the interpolated gene expression matrix and 
#' \code{time.series} the time of the interpGE matrix columns.
#' 
#' @seealso \code{\link[RAPToR]{prepare_refdata}} \code{\link[RAPToR]{ge_im}}
#' 
#' @name Cel_prep
NULL

#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR ge_im
#' @importFrom stats predict
#' 
.prepref_Cel_embryo <- function(n.inter){
  # utils::data("Cel_embryo", envir = environment())
  m <- RAPToR::ge_im(
    X = wormRef::Cel_embryo$g,
    p = wormRef::Cel_embryo$p,
    formula = wormRef::Cel_embryo$geim_params$formula,
    method = wormRef::Cel_embryo$geim_params$method,
    dim_red = wormRef::Cel_embryo$geim_params$dim_red,
    nc = wormRef::Cel_embryo$geim_params$nc
  )
  ndat <- data.frame(age = seq(min(wormRef::Cel_embryo$p$age),
                               max(wormRef::Cel_embryo$p$age),
                               l = n.inter))
  return(
    list(interpGE = predict(m, ndat), time.series = ndat$age)
  )
}

#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR ge_im
#' @importFrom stats predict
#'
.prepref_Cel_larval <- function(n.inter){
  # utils::data("Cel_larval", envir = environment())
  m <- RAPToR::ge_im(
    X = wormRef::Cel_larval$g,
    p = wormRef::Cel_larval$p,
    formula = wormRef::Cel_larval$geim_params$formula,
    method = wormRef::Cel_larval$geim_params$method,
    dim_red = wormRef::Cel_larval$geim_params$dim_red,
    nc = wormRef::Cel_larval$geim_params$nc
  )
  ndat <- data.frame(age = seq(min(wormRef::Cel_larval$p$age),
                               max(wormRef::Cel_larval$p$age),
                               l = n.inter),
                     cov = rep(wormRef::Cel_larval$p$cov[1], n.inter))
  return(
    list(interpGE = predict(m, ndat), time.series = ndat$age)
  )
}


#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR ge_im
#' @importFrom stats predict
#' 
.prepref_Cel_larv_YA <- function(n.inter){
  # utils::data("Cel_larv_YA", envir = environment())
  m <- RAPToR::ge_im(
    X = wormRef::Cel_larv_YA$g,
    p = wormRef::Cel_larv_YA$p,
    formula = wormRef::Cel_larv_YA$geim_params$formula,
    method = wormRef::Cel_larv_YA$geim_params$method,
    dim_red = wormRef::Cel_larv_YA$geim_params$dim_red,
    nc = wormRef::Cel_larv_YA$geim_params$nc
  )
  ndat <- data.frame(age = seq(min(wormRef::Cel_larv_YA$p$age),
                               max(wormRef::Cel_larv_YA$p$age),
                               l = n.inter))
  return(
    list(interpGE = predict(m, ndat), time.series = ndat$age)
  )
}

#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR ge_im
#' @importFrom stats predict
#'
.prepref_Cel_YA_1 <- function(n.inter){
  # utils::data("Cel_YA_1", envir = environment())
  m <- RAPToR::ge_im(
    X = wormRef::Cel_YA_1$g,
    p = wormRef::Cel_YA_1$p,
    formula = wormRef::Cel_YA_1$geim_params$formula,
    method = wormRef::Cel_YA_1$geim_params$method,
    dim_red = wormRef::Cel_YA_1$geim_params$dim_red,
    nc = wormRef::Cel_YA_1$geim_params$nc
  )
  ndat <- data.frame(age = seq(min(wormRef::Cel_YA_1$p$age),
                               max(wormRef::Cel_YA_1$p$age),
                               l = n.inter),
                     cov = rep("N2.NI", n.inter))
  return(
    list(interpGE = predict(m, ndat), time.series = ndat$age)
  )
}

#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR ge_im
#' @importFrom utils data
#'
.prepref_Cel_YA_2 <- function(n.inter){
  # utils::data("Cel_YA_2", envir = environment())
  m <- RAPToR::ge_im(
    X = wormRef::Cel_YA_2$g,
    p = wormRef::Cel_YA_2$p,
    formula = wormRef::Cel_YA_2$geim_params$formula,
    method = wormRef::Cel_YA_2$geim_params$method,
    dim_red = wormRef::Cel_YA_2$geim_params$dim_red,
    nc = wormRef::Cel_YA_2$geim_params$nc
  )
  ndat <- data.frame(age = seq(48, # bc of bad quality
                               max(wormRef::Cel_YA_2$p$age),
                               l = n.inter),
                     cov = rep(wormRef::Cel_YA_2$p$cov[2], n.inter))
  return(
    list(interpGE = predict(m, ndat), time.series = ndat$age)
  )
}