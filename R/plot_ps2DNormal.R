#' Plotting function for \code{ps2DNormal}
#'
#' @description Plotting function for 2D P-spline smooothing
#' (using \code{ps2DNormal} with \code{class ps2dnormal}).
#' @import graphics
#' @param x the P-spline object, usually from \code{ps2DNormal}.
#' @param ... other parameters.
#' @param xlab label for the x-axis, e.g. "my x" (quotes required).
#' @param ylab label for the y-axis, e.g. "my y" (quotes required).
#' @param Resol resolution for plotting, default \code{Resol = 100}.

#' @return
#' \item{Plot}{a plot of the smooth 2D P-spline smooth surface.}
#'
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#' @examples
#' library(SemiPar)
#' library(fields)
#' library(spam)
#' library(JOPS)
#'
#' # Get the data
#' data(ethanol)
#' x <- ethanol$C
#' y <- ethanol$E
#' z <- ethanol$NOx
#'
#' # Set parameters for domain
#' xlo <- 7
#' xhi <- 19
#' ylo <- 0.5
#' yhi <- 1.25
#'
#' # Set P-spline parameters, fit and compute surface
#' xpars <- c(xlo, xhi, 10, 3, 3, 1)
#' ypars <- c(ylo, yhi, 10, 3, 3, 1)
#' Pars1 <- rbind(xpars, ypars)
#' fit <- ps2DNormal(cbind(x, y, z), Pars = Pars1)
#' plot(fit, xlab = "C", ylab = "E")
#' @export

#library(fields)
plot.ps2dnormal <- function(x,..., xlab = " ", ylab = " ", Resol = 100) {
  ps2dnor = x
  Pars <- ps2dnor$Pars
  Data <- ps2dnor$Data
  x <- Data[, 1]
  y <- Data[, 2]
  z <- Data[, 3]

  Xgrid <- seq(Pars[1, 1], Pars[1, 2], length = Resol)
  Ygrid <- seq(Pars[2, 1], Pars[2, 2], length = Resol)
  XYgrid <- as.matrix(expand.grid(Xgrid, Ygrid))
  fit <- ps2DNormal(cbind(x, y, z), Pars, XYpred = XYgrid)
  image.plot(Xgrid, Ygrid, matrix(fit$pred, Resol, Resol), xlab = xlab, ylab = ylab)
}
