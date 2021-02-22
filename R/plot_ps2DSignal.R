#' Plotting function for \code{ps2DSignal}
#'
#' @description Plotting function for 2D P-spline signal regression
#' coefficients (using \code{ps2DSignal} with \code{class ps2dsignal}). Although
#' standard error surface bands
#' can be comuputed they are intentially left out as they are not
#' interpretable, and there is generally little data to steer
#' such a high-dimensional parameterization.
#'
#' @import fields
#'
#' @param x the P-spline object, usually from \code{ps2DSignal}.
#' @param ... other parameters.
#' @param xlab label for the x-axis, e.g. "my x" (quotes required).
#' @param ylab label for the y-axis, e.g. "my y" (quotes required).
#' @param Resol Resolution of bgrid (default \code{Resol = 200}).
#' @return
#' \item{Plot}{a plot of the 2D P-spline signal coefficent surface.}
#'
#' @author Paul Eilers and Brian Marx
#' @references Marx, B.D. and Eilers, P.H.C. (2005).
#' Multidimensional penalized signal regression, \emph{Technometrics}, 47: 13-22.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' library(fields)
#' library(JOPS)
#'
#' # Get the data
#' x0 <- Sugar$X
#' x0 <- x0 - apply(x0, 1, mean) # center Signal
#' y <- as.vector(Sugar$y[, 3]) # Response is Ash
#'
#' # Inputs for two-dimensional signal regression
#' nseg <- c(7, 37)
#' pord <- c(3, 3)
#' min_ <- c(230, 275)
#' max_ <- c(340, 560)
#' M1_index <- rev(c(340, 325, 305, 290, 255, 240, 230))
#' M2_index <- seq(from = 275, to = 560, by = .5)
#' p1 <- length(M1_index)
#' p2 <- length(M2_index)
#'
#' # Fit optimal model based on LOOCV
#' opt_lam <- c(8858.6679, 428.1332) # Found via svcm
#' Pars_opt <- rbind(
#'  c(min_[1], max_[1], nseg[1], 3, opt_lam[1], pord[1]),
#'  c(min_[2], max_[2], nseg[2], 3, opt_lam[2], pord[2]))
#'
#' fit <- ps2DSignal(y, x0, p1, p2, "unfolded", M1_index, M2_index,
#'        Pars_opt, int = FALSE, ridge_adj = 1e-4 )
#'
#' # Plotting coefficient image
#' plot(fit)
#' @export
#'
# library(fields)
plot.ps2dsignal <- function(x,..., xlab = " ", ylab = " ",
                            Resol = 200) {

  # Prepare estimated alpha surface
  ps2dsig = x
  M1index <- ps2dsig$M1index
  M2index <- ps2dsig$M2index
  pcoef <- ps2dsig$pcoef
  Pars <- ps2dsig$Pars
  int = ps2dsig$int

  # High resolution Grid Bases
  M1grid <- seq(from = min(M1index), to = max(M1index), length = Resol)
  M2grid <- seq(from = min(M2index), to = max(M2index), length = Resol)
  p1 <- length(M1grid)
  p2 <- length(M2grid)
  oM1g <- outer(rep(1, p2), M1grid)
  B1grid <- bbase(as.vector(oM1g), Pars[1, 1], Pars[1, 2], Pars[1, 3], Pars[1, 4])
  oM2g <- outer(M2grid, rep(1, p1))
  B2grid <- bbase(as.vector(oM2g), Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
  n1g <- ncol(B1grid)
  n2g <- ncol(B2grid)

  # Compute tensor products for estimated alpha surface
  B1_g <- kronecker(B1grid, t(rep(1, n2g)))
  B2_g <- kronecker(t(rep(1, n1g)), B2grid)
  Bgrid <- B1_g * B2_g
  ncolB <- ncol(Bgrid)
  A_hat_grid <- Bgrid %*% pcoef[(1 + int):(length(pcoef))]
  A_hatm_grid <- matrix(A_hat_grid, Resol, Resol, byrow = TRUE)


  # Image plot of smooth coefficents
  image.plot(M2grid, M1grid, t(A_hatm_grid),
    xlab = xlab, ylab = ylab, sub = "Smooth Coefficient Surface"
  )

}
