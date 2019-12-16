#' Build Interpolated Gene Expression Reference
#' 
#' Builds the PLSR interpolation of the reference dataset.
#' **These functions are internally called by \code{\link[RAPToR]{prepare_refdata}} from RAPToR.**
#' 
#' @param n.inter passed on to \code{\link[RAPToR]{plsr_interpol}}
#' 
#' @return The output of \code{\link[RAPToR]{plsr_interpol}}
#' @seealso \code{\link[RAPToR]{plsr_interpol}}
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
      n.inter = n.inter)
    )
}

#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR plsr_interpol
#' @importFrom utils data
#'
.prepref_Cel_larval <- function(n.inter){
  # utils::data("Cel_larval", envir = environment())
  return(
    RAPToR::plsr_interpol(
      X = wormRef::Cel_larval$g, 
      time.series = wormRef::Cel_larval$p$age, 
      covar = wormRef::Cel_larval$p$cov, 
      df = wormRef::Cel_larval$df,
      n.inter = n.inter)
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
      n.inter = n.inter,
      tmin = 48) # because of bad quality
  )
}