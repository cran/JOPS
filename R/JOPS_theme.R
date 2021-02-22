#' Custom theming function used to unify ggplot features
#'
#' @description Set a ggplot theme in black and white, with centered titles.
#' @param h_just horizontal justification for ggplot2.
#' @return Custom theming function used to unify ggplot features.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @export
#' @import ggplot2
JOPS_theme = function(h_just = 0.5) {
  ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = h_just))
}


