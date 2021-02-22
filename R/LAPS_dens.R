#' Bayesian density estimation
#'
#' @description Bayesian density estimation with P-splines and Laplace approximation.
#'
#' @param y vector (length \code{m}) of counts, usually a histogram.
#' @param B matrix (\code{m} by \code{n}) with B-spline basis, see \code{bbase()}.
#' @param P penalty matrix (\code{n} by \code{n}).
#' @param loglambdas a vector of values of logarithms of \code{lambda} to explore.
#' @param tol convergence tolerance (relative change in coefficients), default \code{1e-5}.
#' @param mon TRUE or FALSE to monitor the iteration history (default FALSE).
#'
#' @return A list with elements:
#' \item{alpha}{P-spline coefficients of length \code{n}.}
#' \item{weights}{weights from the Laplace approximation, which sum to 1 and are
#' the same length as \code{loglambdas}.}
#' \item{mu}{a vector of length \code{m} of expected values.}
#' \item{Cov}{covariance matrix (\code{m} by \code{m}) of \code{log(mu)}.}
#' \item{lambda}{the penalty parameter.}
#' \item{ed}{the effective model dimension.}
#'
#' @details
#' The B-spline basis should be based on the midpoints of the histogram bins.
#' See the example below.
#' This function is based on the paper of Gressani and Lambert (2018) and code input by Oswaldo Gressani.
#'
#' @author Paul Eilers
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references Gressani, O. and Lambert, P. (2018).
#' Fast Bayesian inference using Laplace approximations in a flexible promotion time cure model based on P-splines.
#' \emph{Computational Statistics and Data Analysis} 124, 151-167.
#'
#' @examples
#' # Smoothing a histogram of Old Faithful eruption durations
#' data(faithful)
#' durations = faithful[, 1]  # Eruption length
#'
#' # Histogram with narrow bin widths
#' bw = 0.05
#' hst = hist(durations, breaks = seq(1, 6, by = bw), plot = TRUE)
#' x = hst$mids
#' y = hst$counts
#'
#' # B-spline basis matrices, for fitting and plotting
#' nseg = 30
#' B = bbase(x, nseg = nseg)
#' xg = seq(min(x), max(x), by = 0.01)
#' Bg = bbase(xg, nseg = nseg)
#' n = ncol(B)
#'
#' # Penalty matrix
#' D2 = diff(diag(n), diff = 2)
#' P2 = t(D2) %*% D2
#'
#' # Fit the model
#' loglambs = seq(-1, 2, by = 0.05)
#' laps2 = LAPS_dens(B, P2, y, loglambs, mon = FALSE)
#' fhat2 = exp(Bg %*% laps2$alpha)
#' lines(xg, fhat2, col = "blue", lwd = 2)

#' @export

LAPS_dens = function(B, P, y, loglambdas, tol = 1e-5, mon = FALSE) {

  # Prior parameters
  a <- 0.1
  b <- 0.00001

  # Log-likelihood
  loglik <- function(B, theta, y) {
    eta <- B %*% theta
    return(sum(y * eta - exp(eta)))
  }

  # Prepare innner products
  n = ncol(B)
  BtB <- t(B) %*% B
  Bty <- t(B) %*% y
  Thetas <- Vars <- Covs <- lprobs <- eds <- NULL

  # Initalize
  mu <- y + 1
  eta <- log(mu)
  lambdas <- 10 ^ loglambdas
  theta <- solve(t(B) %*% diag(c(mu)) %*% B + lambdas[1] * P, t(B) %*% (mu * eta))
  nl <- length(lambdas)

  for (lambda in lambdas) {
    # Compute posterior mode
    for (it in 1:20) {
      eta <- B %*% theta
      mu <- c(exp(eta))
      W <- diag(mu)
      H <- -t(B) %*% W %*% B
      Prec <- lambda * P - H
      g <- t(B) %*% (y - mu) - H %*% theta
      thetnew <- solve(Prec, g)
      dthet <- max(abs(thetnew - theta))
      theta <- thetnew
      if (mon) cat(lambda, it, dthet, "\n")
      if (dthet < tol) break
    }
    C <- chol(Prec)
    V <- chol2inv(C)
    vars <- diag(V)
    ed = sum(diag(V %*% H))

    # Save results
    Thetas <- cbind(Thetas, theta)
    Vars <- cbind(Vars, vars)
    Covs <- cbind(Covs, as.vector(V))
    eds = c(eds, ed)

    # Compute logs of posterior probabilities
    logclamb <- loglik(B, theta, y) - sum(log(diag(C)))
    u2 <- -lambda * t(theta) %*% P %*% theta / 2 + b
    logu1 <- (n / 2 + a - 1) * log(lambda)
    loggammval <- logu1 + u2
    logprob <- loggammval + logclamb
    lprobs <- c(lprobs, logprob)
  }

  # Compute weights
  lprobs <- lprobs - max(lprobs)
  weights <- exp(lprobs)
  weights <- weights / sum(weights)

  # Compute weighted averages
  lla_mean = sum(loglambdas * weights)
  lambda_mean = exp(lla_mean)
  ed_mean = sum(eds * weights)
  theta_mean <- Thetas %*% weights
  var_mean <- Vars %*% weights
  # lambda_mean = lambdas %*% weights
  lambda_mean <- exp(log(lambdas) %*% weights)
  Cov_mean <- matrix(Covs %*% weights, n, n)
  Thc = Thetas - outer(c(theta_mean), rep(1, nl))
  for (k in 1:nl) Cov_mean = Cov_mean + weights[k] * outer(Thc[, k], Thc[, k] )
  sdtheta = sqrt(diag(Cov_mean))

  res = list(alpha = theta_mean, weights = weights, Cov = Cov_mean,
             mu = mu, lambda = lambda_mean, ed = ed_mean)

}
