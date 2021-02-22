#' Smoothing scattered (normal) data using P-splines.
#'
#' @description \code{psNormal} is used to smooth scattered (normal) data using P-splines (with identity link function).
#'
#' @import stats
#' @param y the response vector, usually continuous data.
#' @param x the vector for the continuous regressor of \code{length(y)} and the abcissae used to build the B-spline basis.
#' @param wts the vector of general weights, default is 1; zero allowed.
#' @param xl the number for the min along \code{x} (default is min(\code{x})) .
#' @param xr the number for the max along \code{x} (default is max(\code{x})).
#' @param nseg the number of evenly spaced segments between \code{xl} and \code{xr}.
#' @param bdeg the number of the degree of the basis, usually 1, 2 (default), or 3.
#' @param pord the number of the order of the difference penalty, usually 1, 2, or 3 (defalult).
#' @param lambda the (positive) number for the tuning parameter for the penalty (default 1).
#' @param xgrid a scalar or a vector that gives the \code{x} locations for prediction, useful for plotting.
#' If a scalar (default 100) is used then a uniform grid of this size along (\code{xl}, \code{xr}).

#'
#' @return
#' \item{pcoeff}{a vector of length \code{n} of estimated P-spline coefficients.}
#' \item{muhat}{a vector of length \code{m} of smooth estimated means.}
#' \item{B}{a matrix of dimension \code{m} by \code{n} for the B-spline basis matrix.}
#' \item{wts}{a vector of length \code{m} of weights.}
#' \item{effdim}{estimated effective dimension.}
#' \item{ed_resid}{approximate df residual.}
#' \item{sigma}{square root of MSE.}
#' \item{cv}{standard error of leave-one-out prediction or root average PRESS.}
#' \item{nseg}{the number of B-spline segments.}
#' \item{bdeg}{the degree of the B-spline basis.}
#' \item{pord}{the order of the difference penalty.}
#' \item{lambda}{the positive tuning parameter.}
#' \item{xgrid}{gridded x values, useful for plotting.}
#' \item{ygrid}{gridded fitted mean values, useful for plotting.}
#' \item{se_eta}{gridded standard errors for the fitted mean values, useful for plotting.}
#' \item{P}{"half" of the penalty, such that \code{P'P= lambda D'D}.}
#'
#' @author  Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references  Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#' @examples
#' library(JOPS)
#' library(MASS)
#' data(mcycle)
#' x <- mcycle$times
#' y <- mcycle$accel
#' fit1 <- psNormal(x, y, nseg = 20, bdeg = 3, pord = 2, lambda = .8)
#' plot(fit1, se = 2, xlab = "Time (ms)", ylab = "Acceleration")
#' @export
#'
psNormal <- function(x, y, xl = min(x), xr = max(x), nseg = 10, bdeg = 3,
                     pord = 2, lambda = 1, wts = NULL, xgrid = 100) {
  m <- length(x)
  B <- bbase(x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)

  # Construct penalty stuff
  n <- dim(B)[2]
  P <- sqrt(lambda) * diff(diag(n), diff = pord)
  nix <- rep(0, n - pord)

  # Fit
  if (missing(wts)) {
    wts <- rep(1, m)
  }
  f <- lsfit(rbind(B, P), c(y, nix), intercept = FALSE, wt = c(wts, (nix + 1)))
  qr <- f$qr
  h <- hat(qr)[1:m]
  beta <- f$coef
  mu <- B %*% beta

  # Cross-validation and dispersion
  r <- (y - mu) / (1 - h)
  cv <- sqrt(mean(r^2))
  ed <- sum(h)
  sigma <- sqrt(sum((y - mu)^2) / (m - ed))

  # Compute curve on grid
  if (length(xgrid) == 1) {
    xgrid <- seq(xl, xr, length = xgrid)
  }
  Bu <- bbase(xgrid, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)
  zu <- Bu %*% beta
  ygrid <- zu

  # SE bands on a grid using QR
  R <- qr.R(qr)
  L <- forwardsolve(t(R), t(Bu))
  v2 <- sigma^2 * colSums(L * L)
  se_eta <- sqrt(v2)

  # Return list
  pp <- list(
    x = x, y = y, B = B, P = P, muhat = mu, nseg = nseg, xl = xl,
    xr = xr, bdeg = bdeg, pord = pord, lambda = lambda,
    cv = cv, effdim = ed, ed_resid = m - ed, wts = wts,
    pcoeff = beta, family = "gaussian", link = "identity",
    sigma = sigma, xgrid = xgrid, ygrid = ygrid, se_eta = se_eta
  )
  class(pp) <- "pspfit"
  return(pp)
}
