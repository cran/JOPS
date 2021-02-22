#' Translated number vector to bin index.
#'
#' @description Translates number vector to bin index, given lower and
#' upper limits of the domain and number of bins. A support function
#' for (smoothing) histograms.
#'
#' @param x a numerical vector.
#' @param xmin the lower limit of the domain.
#' @param xmax the upper limit of the domain.
#' @param nbin the number of bins (default=100).
#' @return A list with components:
#' \item{xbin}{a vector of \code{length(x)} with elements giving the bin index.}
#' \item{xgrid}{a vector of \code{length(nbin)} with the midpoints of the bins.}
#' \item{nbin}{the number of bins.}
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @export

binit = function(x, xmin = min(x), xmax = max(x), nbin = 100) {
  dx <- 1.000001 * (xmax - xmin) / nbin  # Little trick to fight effects of rounding
  xbin <- floor(1 + (x - xmin) / dx)
  xgrid <- xmin + (1:nbin - 0.5) * dx
  b = list(xbin = xbin, xgrid = xgrid, nbin = nbin)
  b }
