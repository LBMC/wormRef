#' Plot Reference Timelines
#' 
#' Plots the coverage of the `wormRef` references (requires ggplot2)
#' 
#' @examples 
#' \donttest{
#' plot_refs_Cel()
#' }
#' 
#' @export
#' 
#' @importFrom utils data
plot_refs_Cel <- function(){
  if(!requireNamespace("ggplot2")){
    message("Install ggplot2 to see this graph !")
  }
  utils::data("Cel_devstages", envir = environment())
  dat <- wormRef::Cel_devstages
  lbs <- data.frame(tunit=unique(dat$devstages$tunit), lab = c("Embryo", "Post-Hatching"))
  
  eb <- ggplot2::element_blank()
    
  ggplot2::ggplot(data = dat$datasets, mapping = ggplot2::aes(x = tstart, y = name, xend = tend,)) +
    ggplot2::geom_segment(data = dat$datasets, ggplot2::aes( yend = name, col = name, size = 1.5), show.legend = F) +
    ggplot2::geom_text(data = dat$datasets, ggplot2::aes(col = name, label = name, size = 2), 
                       fontface = "bold", vjust = "bottom", hjust = "left", nudge_y = .1,
                       show.legend = F) +
    ggplot2::geom_text(data = lbs, ggplot2::aes(label = tunit), 
                       x = Inf, y = -Inf, hjust = 1.5, vjust = -1, inherit.aes = F) +
    
    ggplot2::geom_segment(data = dat$devstages, ggplot2::aes(size = 1, y = "DevStage", yend = "DevStage"), 
                          color = rep(c("grey50", "grey90"), length.out = nrow(dat$devstages)), show.legend = F) +
    ggplot2::geom_text(data = dat$devstages, ggplot2::aes(y = "DevStage", label = name), size = 3, 
                       fontface = "bold", vjust = "bottom", hjust = "left", nudge_y = .1,
                       show.legend = F) +
    
    ggplot2::facet_wrap(~tunit, scales = "free", nrow = 2) +
    ggplot2::theme_classic() + ggplot2::xlab("time") + 
    ggplot2::theme(strip.background = eb, strip.text = eb, 
                   axis.title.y = eb, axis.text.y = eb) 
}
