#' Compute a truncated power function.
#'
#' @param x a vector on which the basis is calculated.
#' @param knot a scalar giving the truncation point.
#' @param p a scalar power for the basis, e.g. \code{p = 3} for cubic TPF.
#'
#' @return a vector with the truncated power function.
#'
#' @author Paul Eilers
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' library(JOPS)
#' # Basis  on grid
#' x = seq(0, 4, length = 500)
#' knots = 0:3
#' Y = outer(x, knots, tpower, 1)
#' matplot(x, Y, type ='l', lwd = 2, xlab = 'x', ylab = '',
#' main ='Linear TPF basis')
#'
#' @export

tpower <- function(x, knot, p) {
  (x - knot) ^ p * (x >= knot)
}

