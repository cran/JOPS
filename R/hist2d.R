#' Compute a 2D histogram
#'
#' @description Compute a two-dimesnional histogram from two vectors (of the same length), \code{x} and \code{y}.
#'
#' @param x a numeric vector.
#' @param y a numeric vector of the same length as \code{x}.
#' @param nb a vector \code{c(nbx, nby)}, or a scalar \code{nb}, providing the number of bins for x, and y; default is 100; see details.
#' @param xlim a vector \code{c(xmin, xmax)} containing the limits of the domain of \code{x}; default \code{range(x)}.
#' @param ylim a vector \code{c(ymin, ymax)} containing the limits of the domain of \code{y}; default \code{range(y)}.
#'
#' @return A list with components:
#' \item{H}{a matrix of dimension \code{nbx} by \code{nby} containing bin counts.}
#' \item{xgrid}{a vector of length \code{nbx} representing centers of the bins for \code{x}.}
#' \item{ygrid}{a vector of length \code{nby} representing centers of the bins for \code{y}.}
#' \item{xbin}{a vector giving the bin number of each element of \code{x}.}
#' \item{ybin}{a vector giving the bin number of each element of \code{y}.}
#'
#' @details
#' If \code{nb} is scalar, it is extended to \code{c(nb, nb)}, so that both dimensions will have the same number of bins.
#'
#' Elements of \code{x} (\code{y}) that fall outside the range specified by \code{xlim} (\code{ylim}) are not counted.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' data(faithful)
#' x = faithful$eruptions
#' y = faithful$waiting
#' C = hist2d(x, y, c(50,50))
#' image(C$xgrid, C$ygrid, C$H, xlab='Eruption length (min)', ylab='Waiting time (min)')
#' title('Old Faithful geyser')
#'
#' @export


hist2d = function(x, y, nb = c(100, 100), xlim = range(x), ylim = range(y)) {
  if (length(nb) == 1) nb = c(nb,nb )
  xb = binit(x, xmin = xlim[1], xmax = xlim[2], nbin = nb[1])
  yb = binit(y, xmin = ylim[1], xmax = ylim[2], nbin = nb[2])
  xx = xb$xbin
  yy = yb$xbin
  sel = 0 < xx & xx <= nb[1] & 0 < yy & yy <= nb[2]
  H = count2d(xx[sel], yy[sel], nb)
  h2d = list(H = H, xgrid = xb$xgrid, ygrid = yb$xgrid, xbin = xx, ybin = yy)
}
