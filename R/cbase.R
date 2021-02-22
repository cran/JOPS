#' Compute a circular B-spline basis matrix
#'
#' @description Computes a circular B-spline basis matrix using evenly spaced knots.
#'
#' @param x a vector of argument values, at which the B-spline basis functions
#' are to be evaluated.
#' @param xl the lower limit of the domain of x; default is \code{min(x)}.
#' @param xr the upper limit of the domain of x; default is \code{max(x)}.
#' @param nseg the number of B-spline segments (default 10) between xl and xr.
#' @param bdeg the degree of the basis, usually 1, 2, or 3 (default).
#'
#' @return A matrix with \code{length(x)} rows and \code{nseg} columns.
#'
#' @details
#' If \code{xl} is larger than \code{min(x)}, it wil be adjusted to \code{min(x)} and a warning wil be given.
#' If \code{xr} is smaller than \code{max(x)}, it wil be adjusted to \code{max(x)} and a warning wil be given.
#'
#' The design parameters \code{x, xl, xr, ndeg, bdeg} and \code{type = 'cbase'} are added to the list of attributes.
#'
#' In a circular basis, the B-splines are wrapped around the boundaries of the domain. Use a circular basis for data
#' like directions or angles. It should be combined with a circular penalty matrix, as computed by \code{cdiff()}.
#'
#' @author Paul Eilers and Brian Marx
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M.
#' (2015). Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#' @examples
#' # Compute and plot a circular B-spline basis matrix
#' x = seq(0, 360, by = 2)
#' B = cbase(x, 0, 360, nseg = 8, bdeg = 3)
#' matplot(x, B, type = 'l', lty = 1, lwd = 2, xlab = 'x', ylab = '')
#' title('Note how the ends connect smoothly meet at boundaries' )
#'
#' @export

cbase <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) {

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
  B0 = bbase(x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)
  n = ncol(B0) - bdeg
  cc = (1:bdeg) + n
  B = B0[, 1:n]
  B[, 1:bdeg] = B[, 1:bdeg] + B0[, cc]

  # Add design information as attributes
  att1 = attributes(B)
  att2 = list(x = x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg, type = 'cbase')
  attributes(B) <- append(att1, att2)
  return(B)
}

