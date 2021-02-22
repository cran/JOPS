#' Themeing functions used to unify ggplot features
#'
#' @description  Custom size and color of points.
#' @param s_size point size parameter for ggplot2 (default = 1.5).
#' @return themeing function for ggplot2 features.
#' @export
#' @import ggplot2

JOPS_point = function(s_size = 1.5) {
  ggplot2::geom_point(colour = grey(0.5), size = s_size)
}

#' Custom theme for ggplot
#'
#' @description Set a ggplot theme in black and white, with centered titles.
#' @param h_just hjust parameter setting for ggplot.
#' @return custom theme for ggplot.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @export
#' @import ggplot2
JOPS_theme = function(h_just = 0.5) {
  ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = h_just))
}

#' Custom color ramp.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @param n number of steps.
#' @return custom color ramp.
#' @export
JOPS_colors = function(n) colorspace::rainbow_hcl(n, start = 10, end = 350)

