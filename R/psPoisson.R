#' Smoothing scattered Poisson data using P-splines.
#'
#' @description \code{psPoisson} is used to smooth scattered
#' Poisson data using P-splines with a log link function.
#'
#' @import stats
#'
#' @param y the response vector, usually count data.
#' @param x the vector for the continuous regressor of \code{length(y)} and
#' the abcissae used to build the B-spline basis.
#' @param xl the number for the min along \code{x} (default is min(\code{x})).
#' @param xr the number for the max along \code{x} (default is max(\code{x})).
#' @param nseg the number of evenly spaced segments between \code{xl} and \code{xr} (default 10).
#' @param bdeg the number of the degree of the basis, usually 1, 2, or 3 (defalult).
#' @param pord the number of the order of the difference penalty, usually 1, 2 (default), or 3.
#' @param lambda the (positive) number for the tuning parameter for the penalty (default 1).
#' @param wts the vector of general weights, zeros are allowed (default 1).
#' @param show Set to TRUE or FALSE to display iteration history (default FALSE).
#' @param iter a scalar to set the maximum number of iterations, default \code{iter=100}.
#' @param xgrid a scalar or a vector that gives the \code{x} locations for prediction, useful for plotting.
#' If a scalar (default 100) is used then a uniform grid of this size along (\code{xl}, \code{xr}).

#'
#' @return
#' \item{pcoef}{a vector of length \code{n} of estimated P-spline coefficients.}
#' \item{muhat}{a vector of length \code{m} of estimated means.}
#' \item{B}{the \code{m} by \code{n} B-spline basis.}
#' \item{dev}{deviance of fit.}
#' \item{effdim}{effective dimension of fit.}
#' \item{aic}{AIC.}
#' \item{wts}{the vector of given prior weights.}
#' \item{nseg}{the number of B-spline segments.}
#' \item{bdeg}{the degree of the B-spline basis.}
#' \item{pord}{the order of the difference penalty.}
#' \item{lambda}{the positive tuning parameter.}
#' \item{family}{the family of the response (
#' \code{"Poisson"}).}
#' \item{link}{the link function used (\code{"log"}).}
#' \item{xgrid}{gridded x values, useful for plotting.}
#' \item{ygrid}{gridded fitted linear predictor values, useful for plotting.}
#' \item{mugrid}{gridded (inverse link) fitted mean values, useful for plotting.}
#' \item{se_eta}{gridded standard errors for the linear predictor.}
#' \item{dispersion}{Dispersion parameter estimated \code{dev/(m-effdim)}.}
#'
#'
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#' @examples
#' library(JOPS)
#' library(boot)
#'
#' # Extract the data
#' Count <- hist(coal$date, breaks = c(1851:1963), plot = FALSE)$counts
#' Year <- c(1851:1962)
#' xl <- min(Year)
#' xr <- max(Year)
#'
#' # Poisson smoothing
#' nseg <- 20
#' bdeg <- 3
#' fit1 <- psPoisson(Year, Count, xl, xr, nseg, bdeg, pord = 2, lambda = 1)
#' plot(fit1, xlab = "Year", ylab = "Count", se = 2)
#' @export

psPoisson <- function(x, y, xl = min(x), xr = max(x), nseg = 10, bdeg = 3,
                      pord = 2, lambda = 1, wts = NULL, show = FALSE, iter = 100, xgrid = 100) {

  # Compute B-spline basis
  m <- length(x)
  B <- bbase(x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)

  # Construct penalty stuff
  n <- dim(B)[2]
  P <- sqrt(lambda) * diff(diag(n), diff = pord)
  nix <- rep(0, n - pord)

  # Initialize
  z <- log(y + 0.01)
  if (missing(wts)) {
    wts <- rep(1, m)
  }
  # Fit
  for (it in 1:iter) {
    mu <- exp(z)
    w <- mu
    u <- (y - mu) / w + z
    wtprod <- c(wts * w, (nix + 1))
    f <- lsfit(rbind(B, P), c(u, nix), intercept = FALSE, wt = wtprod)
    beta <- f$coef
    znew <- B %*% beta
    dz <- max(abs(z - znew))
    z <- znew
    if (dz < 1e-06) {
      break
    }
    if (show) {
      print(c(it, dz))
    }
  }

  if (it > (iter - 1)) {
    cat(paste("Did NOT converge, iter >", iter))
    warning(paste("Did NOT converge"))
  }

  # Compute AIC, ED, Dispersion
  dev <- 2 * sum(y * log((y + 1e-09) / mu))
  qr <- f$qr
  h <- hat(qr)[1:m]
  ed <- sum(h)
  aic <- dev + 2 * ed
  dispersion <- dev / (m - ed)

  # Compute curve on grid
  if (length(xgrid) == 1) {
    xgrid <- seq(xl, xr, length = xgrid)
  }
  Bu <- bbase(xgrid, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)
  zu <- Bu %*% beta
  ygrid <- zu
  mugrid <- exp(zu)

  # SE bands on a grid using QR
  R <- qr.R(qr)
  #C2 <- chol2inv(R)
  L <- forwardsolve(t(R), t(Bu))
  v2 <- colSums(L * L)
  se_eta <- sqrt(v2)

  # Return list
  pp <- list(
    xl = xl, xr = xr, aic = aic, x = x, y = y, B = B, P = P,
    muhat = mu, nseg = nseg, bdeg = bdeg, dev = dev,
    pord = pord, pcoef = beta, lambda = lambda,
    effdim = ed, dispersion = dispersion, aic = aic,
    family = "poisson", link = "log", wts = wts,
    xgrid = xgrid, ygrid = ygrid, mugrid = mugrid, se_eta = se_eta
  )
  class(pp) <- "pspfit"
  return(pp)
}
