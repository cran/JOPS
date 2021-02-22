#' Smoothing scattered binomial data using P-splines.
#' @description \code{psBinomial} is used to smooth scattered
#' binomial data using P-splines using a logit link function.
#'
#' @import stats
#'
#' @param y the response vector, usually 0/1 or binomial counts.
#' @param x the vector for the continuous regressor of \code{length(y)} and
#' the abcissae, on which the B-spline basis is constructed.
#' @param ntrials the vector for the number of binomial trials (default = 1).
#' @param xl the lower limit for the domain of \code{x} (default is min(\code{x})).
#' @param xr the upper limit for the domain of \code{x} (default is max(\code{x})).
#' @param nseg the number of evenly spaced segments between xl and xr.
#' @param bdeg the number of the degree of the basis, usually 1, 2 (default), or 3.
#' @param pord the number of the order of the difference penalty, usually 1, 2, or 3 (defalult).
#' @param lambda the (positive) number for the tuning parameter for the penalty.
#' @param wts the vector of weights, default is 1, zeros allowed.
#' @param show Set to TRUE or FALSE to display iteration history.
#' @param iter a scalar to set the maximum number of iterations, default \code{iter = 100}.
#' @param xgrid a scalar or a vector that gives the \code{x} locations for prediction, useful for plotting.
#' If a scalar (default 100) is used then a uniform grid of this size along (\code{xl}, \code{xr}).


#'
#' @return
#' \item{pcoef}{a vector of length \code{n} of estimated P-spline coefficients.}
#' \item{p}{a vector of length \code{m} of estimated probabilities.}
#' \item{muhat}{a vector of length \code{m} of estimated means (\code{ntrials*p}).}
#' \item{dev}{deviance}
#' \item{effdim}{effective dimension of the smooth.}
#' \item{aic}{AIC}
#' \item{wts}{a vector of preset weights (default = 1).}
#' \item{nseg}{the number of B-spline segments.}
#' \item{bdeg}{the degree of the B-spline basis.}
#' \item{pord}{the order of the difference penalty.}
#' \item{family}{the GLM family (repsonse distribution).}
#' \item{link}{the link function.}
#' \item{y}{the binomial response.}
#' \item{x}{the regressor on which the basis is constructed.}
#' \item{P}{"half" of the penalty matrix, \code{P'P = lambda*D'D}.}
#' \item{B}{the B-spline basis.}
#' \item{lambda}{the positive tuning parameter.}
#' \item{dispersion}{dispersion parameter estimated \code{dev/(m-effdim)}.}
#' \item{xgrid}{gridded \code{x} values,useful for plotting.}
#' \item{ygrid}{gridded fitted linear predictor values, useful for plotting.}
#' \item{pgrid}{gridded (inverse link) fitted probability values, useful for plotting.}
#' \item{se_eta}{gridded standard errors for the linear predictor.}
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.

#' @examples
#' library(JOPS)
#' # Extract data
#' library(rpart)
#' Kyphosis <- kyphosis$Kyphosis
#' Age <- kyphosis$Age
#' y <- 1 * (Kyphosis == "present") # make y 0/1
#' fit1 <- psBinomial(Age, y,
#'   xl = min(Age), xr = max(Age), nseg = 20,
#'   bdeg = 3, pord = 2, lambda = 10
#' )
#' names(fit1)
#' plot(fit1, xlab = "Age", ylab = "0/1", se = 2)
#' @export

psBinomial <- function(x, y, xl = min(x), xr = max(x), nseg = 10, bdeg = 3,
                       pord = 2, lambda = 1, ntrials = 0 * y + 1, wts = NULL, show = FALSE, iter = 100, xgrid = 100) {

  # Compute B-spline basis
  m <- length(x)
  B <- bbase(x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)
  # Construct penalty stuff
  n <- dim(B)[2]
  P <- sqrt(lambda) * diff(diag(n), diff = pord)
  nix <- rep(0, n - pord)

  # Initialize
  znew <- 0.5 # log(y + 0.01)
  z <- rep(1, m)
  if (missing(wts)) {
    wts <- rep(1, m)
  }

  # Fit
  for (it in 1:iter) {
    p <- exp(z) / (1 + exp(z))
    mu <- ntrials * p
    w <- ntrials * p * (1 - p)
    u <- (y - mu) / w + z
    wtprod <- c(wts * w, (nix + 1))
    f <- lsfit(rbind(B, P), c(u, nix), wt = wtprod, intercept = FALSE)
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

  # Compute AIC
  e <- 1e-09
  dev <- 2 * sum((y + e) * log((y + e) / mu) +
    (ntrials - y + e) * log((ntrials - y + e) / (ntrials - mu)))
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
  pgrid <- exp(zu) / (1 + exp(zu))

  # SE bands on a grid using QR
  R <- qr.R(qr)
  #C2 <- chol2inv(R)
  L <- forwardsolve(t(R), t(Bu))
  v2 <- colSums(L * L)
  se_eta <- sqrt(v2)

  # Return list
  pp <- list(
    aic = aic, B = B, P = P, wts = wts, xl = xl, xr = xr, x = x,
    y = y, p = p, muhat = mu, nseg = nseg, bdeg = bdeg,
    pord = pord, pcoef = beta, lambda = lambda,
    effdim = ed, dispersion = dispersion, dev = dev,
    family = "binomial", link = "logit", ntrials = ntrials,
    xgrid = xgrid, ygrid = ygrid, pgrid = pgrid, se_eta = se_eta
  )
  class(pp) <- c("pspfit")
  return(pp)
}
