#' Compute a B-spline basis matrix
#'
#' @description Compute a B-spline basis matrix using evenly spaced knots.
#'
#' @param x a vector of argument values, at which the B-spline basis functions
#' are to be evaluated.
#' @param xl the lower limit of the domain of x; default is \code{min(x)}.
#' @param xr the upper limit of the domain of x; default is \code{max(x)}.
#' @param nseg the number of equally sized segments between xl and xr; default is 10.
#' @param bdeg the degree of the splines, usually 1, 2, or 3 (default).
#'
#' @return A matrix with  \code{length(x)} rows and \code{nseg + bdeg} columns.
#'
#' @details
#' If \code{xl} is larger than \code{min(x)}, it will be adjusted to \code{min(x)} and a warning wil be given.
#' If \code{xr} is smaller than \code{max(x)}, it will be adjusted to \code{max(x)} and a warning wil be given.
#' The values of the design parameters \code{x, xl, xr, ndeg, bdeg} and \code{type = 'bbase'} are added to the list of attributes of the matrix.
#'
#' @author Paul Eilers and Brian Marx
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with
#' B-splines and penalties (with comments and rejoinder), \emph{Statistical Science},
#' 11: 89-121.
#'
#' @references Eilers, P.H.C. and B.D. Marx (2010).
#' Splines, knots and penalties. Wiley Interdisciplinary
#' Reviews: Computational Statistics. Wiley: NY. DOI: 10.1002/wics.125
#'
#' @examples
#' # Compute and plot a B-spline basis matrix
#' x = seq(0, 360, by = 2)
#' B = bbase(x, 0, 360, nseg = 8, bdeg = 3)
#' matplot(x, B, type = 'l', lty = 1, lwd = 2, xlab = 'x', ylab = '')
#'
#' @export

bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) {

  # Check domain and adjust it if necessary
  if (xl > min(x)) {
    xl = min(x)
    warning("Left boundary adjusted to min(x) = ", xl)
  }
  if (xr <  max(x)) {
    xr = max(x)
    warning("Right boundary adjusted to max(x) = ", xr)
  }

  # Compute the B-splines
  dx <- (xr - xl) / nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)

  # Make B-splines exactly zero beyond their end knots
  nb <- ncol(B)
  sk <- knots[(1:nb) + bdeg + 1]
  Mask <- outer(x, sk, '<')
  B <- B * Mask

  # Add design information as attributes
  att1 = attributes(B)
  att2 = list(x = x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg, type = 'bbase')
  attributes(B) <- append(att1, att2)

  return(B)
}

