#' Plotting function for \code{sim_psr}
#'
#' @description Plotting function for single-index signal
#' regression with tensor product P-splines (using \code{sim_psr} with \code{class simpsr}).
#'
#' @import graphics
#' @import fields
#'
#' @param x the P-spline object, usually from \code{sim_psr}.
#' @param ... other parameters.
#' @param xlab label for the x-axis, e.g. "my x" (quotes required).
#' @param ylab label for the y-axis, e.g. "my y" (quotes required).
#' @param Resol resolution for plotting, default \code{Resol = 100}.
#'
#' @return
#' \item{Plot}{a two panel plot, one for the estimated P-spline signal coefficent vector, and another for
#' the estimated (unkown) P-spline smooth link function.}
#'
#' @author Paul Eilers, Brian Marx, and Bin Li
#' @references Eilers, P.H.C., B. Li, B.D. Marx (2009).
#' Multivariate calibration with single-index signal regression,
#' \emph{Chemometrics and Intellegent Laboratory Systems}, 96(2), 196-202.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' library(JOPS)
#' # Get the data
#' library(fds)
#' data(nirc)
#' iindex <- nirc$x
#' X <- nirc$y
#' sel <- 50:650 # 1200 <= x & x<= 2400
#' X <- X[sel, ]
#' iindex <- iindex[sel]
#' dX <- diff(X)
#' diindex <- iindex[-1]
#' y <- as.vector(labc[1, 1:40])
#' oout <- 23
#' dX <- t(dX[, -oout])
#' y <- y[-oout]
#'
#' pords <- c(2, 2)
#' nsegs <- c(27, 7)
#' bdegs = c(3, 3)
#' lambdas <- c(1e-6, .1)
#' max_iter <- 100
#'
#' # Single-index model
#' fit <- sim_psr(y, dX, diindex, nsegs, bdegs, lambdas, pords,
#'              max_iter)
#' plot(fit, xlab = "Wavelength (nm)", ylab = " ")
#'
#' @export
#'
# library(fields)
plot.simpsr <- function(x,..., xlab = " ", ylab = " ",
                          Resol = 100) {
  sim = x
  eta <- sim$eta
  f = sim$f
  nsegs <- sim$nsegs
  pords <- sim$pords
  lambdas <- sim$lambdas
  x_index = sim$x_index
  y <- sim$y
  bdegs <- sim$bdegs
  alpha <- sim$alpha

  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(2, 1))

  # Prepare bases for estimated alpha surface
  minS = min(x_index)
  maxS = max(x_index)
  S_index <- seq(from = minS, to = maxS, length = Resol)
  B <- bbase(S_index, minS, maxS, nsegs[1], bdegs[1])

  # Compute tensor products for estimated alpha surface

  a_hat <- B %*% alpha
  plot(S_index, a_hat, type='l',
    xlab = xlab, ylab = ylab,
    main = "Coefficient vector" )
  abline(0, 0, lty= 2, lwd = 2)

  mineta = min(eta)
  maxeta = max(eta)
  eta_index <- seq(from = mineta, to = maxeta, length = Resol)
  feta <- psNormal(eta, y, mineta, maxeta, nsegs[2], bdegs[2],
                   pords[2], lambdas[2], xgrid = sort(eta_index))
  plot(eta_index, feta$ygrid, type='l',
        xlab = "Linear predictor (eta)", ylab = "f(eta)", main = "Link function")
  abline(0, 1, lty = 2, col = 1, lwd = 2)

}
