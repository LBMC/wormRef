#' Build Interpolated Gene Expression Reference
#' 
#' Builds the PLSR interpolation of the reference dataset.
#' **These functions are internally called by \code{\link[RAPToR]{prepare_refdata}} from RAPToR.**
#' 
#' @param n.inter passed on to \code{\link[RAPToR]{plsr_interpol}}
#' 
#' @return A list with \code{interpGE} the interpolated gene expression matrix and 
#' \code{time.series} the time of the interpGE matrix columns.
#' 
#' @seealso \code{\link[RAPToR]{plsr_interpol}} \code{\link[RAPToR]{ge_im}}
#' 
#' @name Cel_prep
NULL

#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR plsr_interpol
#' @importFrom utils data
#' 
.prepref_Cel_embryo <- function(n.inter){
  # utils::data("Cel_embryo", envir = environment())
  return(
    RAPToR::plsr_interpol(
      X = wormRef::Cel_embryo$g, 
      time.series = wormRef::Cel_embryo$p$age, 
      covar = wormRef::Cel_embryo$p$cov, 
      df = wormRef::Cel_embryo$df,
      plsr.nc = wormRef::Cel_embryo$nc,
      n.inter = n.inter)
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
#' @importFrom RAPToR plsr_interpol
#' @importFrom utils data
#'
.prepref_Cel_YA_1 <- function(n.inter){
  # utils::data("Cel_YA_1", envir = environment())
  return(
    RAPToR::plsr_interpol(
      X = wormRef::Cel_YA_1$g, 
      time.series = wormRef::Cel_YA_1$p$age, 
      covar = wormRef::Cel_YA_1$p$cov, 
      df = wormRef::Cel_YA_1$df,
      plsr.nc = wormRef::Cel_YA_1$nc,
      n.inter = n.inter)
  )
}

#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR plsr_interpol
#' @importFrom utils data
#'
.prepref_Cel_YA_2 <- function(n.inter){
  # utils::data("Cel_YA_2", envir = environment())
  return(
    RAPToR::plsr_interpol(
      X = wormRef::Cel_YA_2$g, 
      time.series = wormRef::Cel_YA_2$p$age, 
      covar = wormRef::Cel_YA_2$p$cov, 
      df = wormRef::Cel_YA_2$df,
      plsr.nc = wormRef::Cel_YA_2$nc,
      n.inter = n.inter,
      tmin = 48) # because of bad quality
  )
}