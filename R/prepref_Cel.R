#' Build Interpolated Gene Expression References
#' 
#' Builds the interpolation of the reference datasets.
#' **These functions are internally called by \code{\link[RAPToR]{prepare_refdata}} from RAPToR.**
#' 
#' @param n.inter,by.inter the resolution of the interpolation, as in \code{seq(start, end, length.out = n.inter)} or \code{seq(start, end, by = by.inter)}.
#' 
#' @return A list with \code{interpGE} the interpolated gene expression matrix and 
#' \code{time.series} the time of the interpGE matrix columns.
#' 
#' @seealso \code{\link[RAPToR]{prepare_refdata}} \code{\link[RAPToR]{ge_im}} \code{\link[RAPToR]{make_ref}}
#' 
#' @name Cel_prep
NULL



#' @rdname Cel_prep
#' @export
#' @importFrom RAPToR ge_im make_ref
#' @importFrom stats predict
#' 
.prepref_skel <- function(data, from=NULL, to=NULL){
  # .prepref function factory
  f <- function(n.inter=NULL, by.inter=NULL){
    m <- RAPToR::ge_im(
      X = data$g,
      p = data$p,
      formula = data$geim_params$formula,
      method = data$geim_params$method,
      dim_red = data$geim_params$dim_red,
      nc = data$geim_params$nc
    )
    return(RAPToR::make_ref(m, 
                            n.inter = n.inter,
                            by.inter = by.inter,
                            from = from, 
                            to = to,
                            t.unit = data$t.unit,
                            cov.levels = data$cov.levels,
                            metadata = data$metadata)
    )
  }
  return(f)
}



#' @rdname Cel_prep
#' @export
#' 
.prepref_Cel_embryo <- .prepref_skel(wormRef::Cel_embryo)


#' @rdname Cel_prep
#' @export
#'
.prepref_Cel_larval <- .prepref_skel(wormRef::Cel_larval)


#' @rdname Cel_prep
#' @export
#' 
.prepref_Cel_larv_YA <- .prepref_skel(wormRef::Cel_larv_YA)


#' @rdname Cel_prep
#' @export
#'
.prepref_Cel_YA_1 <- .prepref_skel(wormRef::Cel_YA_1)


#' @rdname Cel_prep
#' @export
#'
.prepref_Cel_YA_2 <- .prepref_skel(wormRef::Cel_YA_2, from=48) 
# from = 48h because of bad quality


#' @rdname Cel_prep
#' @export
#'
.prepref_Cel_aging_1 <- .prepref_skel(wormRef::Cel_aging_1) 
# from = 48h because of bad quality
