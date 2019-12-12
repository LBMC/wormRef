#' Build Interpolated Gene Expression Reference
#' 
#' Builds the PLSR interpolation of the reference dataset.
#' **These functions are internally called by \code{\link[wormAge]{prepare_refdata}} from wormAge.**
#' 
#' @param n.inter passed on to \code{\link[wormAge]{plsr_interpol}}
#' 
#' @return The output of \code{\link[wormAge]{plsr_interpol}}
#' @seealso \code{\link[wormAge]{plsr_interpol}}
#' 
#' @name Cel_prep
NULL

#' @rdname Cel_prep
#' @importFrom wormAge plsr_interpol
#' 
.prepref_Cel_embryo <- function(n.inter){
  return(
    wormAge::plsr_interpol(
      X = Cel_embryo$g, 
      time.series = Cel_embryo$p$age, 
      covar = Cel_embryo$p$cov, 
      df = Cel_embryo$df,
      n.inter = n.inter)
    )
}

#' @rdname Cel_prep
#' @importFrom wormAge plsr_interpol
#'
.prepref_Cel_larval <- function(n.inter){
  return(
    wormAge::plsr_interpol(
      X = Cel_larval$g, 
      time.series = Cel_larval$p$age, 
      covar = Cel_larval$p$cov, 
      df = Cel_larval$df,
      n.inter = n.inter)
  )
}

#' @rdname Cel_prep
#' @importFrom wormAge plsr_interpol
#'
.prepref_Cel_YA_1 <- function(n.inter){
  return(
    wormAge::plsr_interpol(
      X = Cel_YA_1$g, 
      time.series = Cel_YA_1$p$age, 
      covar = Cel_YA_1$p$cov, 
      df = Cel_YA_1$df,
      n.inter = n.inter)
  )
}

#' @rdname Cel_prep
#' @importFrom wormAge plsr_interpol
#'
.prepref_Cel_YA_2 <- function(n.inter){
  return(
    wormAge::plsr_interpol(
      X = Cel_YA_2$g, 
      time.series = Cel_YA_2$p$age, 
      covar = Cel_YA_2$p$cov, 
      df = Cel_YA_2$df,
      n.inter = n.inter,
      tmin = 48) # because of bad quality
  )
}