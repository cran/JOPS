#' Single-Index signal regression using P-splines
#'
#' @description \code{sim_psr} is a single-index
#' signal regression model that estimates both the signal coefficients
#' vector and the unknown link function using P-splines.
#'
#' @import stats
#'
#' @param y a response vector of length \code{m}, usually continuous.
#' @param X The signal regressors with dimension \code{m} by \code{p}.
#' @param x_index an index of length \code{p} for columns of signal matrix;
#' default is simple sequence, \code{c(1: ncol(X))}.
#' @param nsegs a vector of length 2 containing
#' the number of evenly spaced segments between min and max, for each
#' the coefficient vector and the (unknown) link function,
#' resp. (default \code{c(10, 10)}).
#' @param bdegs a vector of length 2 containing
#' the degree of B-splines, for the coefficient vector and
#' the (unknown) link function, resp. (default cubic or \code{c(3, 3)}).
#' @param lambdas a vector of length 2 containing
#' the positive tuning parameters, for each
#' the coefficient vector and the (unknown) link function, resp. (default \code{c(1, 1)}).
#' @param pords a vector of length 2 containing
#' the difference penalty order, for each
#' the coefficient vector and the (unknown) link function, resp. (default\code{c(2, 2)} ).
#' @param max_iter a scalar for the maximum number of iterations (default 100).

#'
#' @return
#' \item{y}{the response vector of length \code{m}.}
#' \item{alpha}{the P-spline coefficient vector of length \code{(nsegs[1]+bdeg[1])}.}
#' \item{iter}{the number of iterations used for the single-index fit.}
#' \item{yint}{the estimated y-intercept for the single-index model.}
#' \item{B}{the B-spline matrix built along the signal index, using \code{nsegs[1]}, used for the coefficient vector.}
#' \item{Q}{the effective regressors from the \code{psVCSignal} portion of the single-index
#' fit with dimension \code{m} by \code{length(alpha)}.}
#' \item{nsegs}{a vector of length 2 containing
#' the number of evenly spaced segments between min and max, for each
#' the coefficient vector and the link function, resp.}
#' \item{bdegs}{a vector of length 2 containing
#' the degree of B-splines, for each
#' the coefficient vector and the link function, resp.}
#' \item{lambdas}{a vector of length 2 containing
#' the positive tuning parameters, for each
#' the coefficient vector and the link function, resp.}
#' \item{pords}{a vector of length 2 containing
#' the difference penalty order, for each
#' the coefficient vector and the link function, resp.}
#' \item{eta}{the estimated linear predictor for the single-index fit.}
#' \item{cv}{the leave-one-out cross-validation statistic or the standard error of prediction for the single-index fit.}
#' \item{delta_alpha}{change measure in signal-coefficent parameters at
#' convervence.}
#' \item{x_index}{the index of length \code{p} for columns of signal matrix.}
#' \item{f_fit}{the \code{psNormal} object, fitting link function f(\code{eta}).}
#' \item{f_eta}{the predicted values of the link function estimated with \code{f_fit} or estimated f(\code{eta}), at \code{x = eta}.}


#' @author Paul Eilers, Brian Marx, and Bin Li
#' @references Eilers, P.H.C., B. Li, B.D. Marx (2009).
#' Multivariate calibration with single-index signal regression,
#' \emph{Chemometrics and Intellegent Laboratory Systems}, 96(2), 196-202.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
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
#'
#' plot(fit, xlab = "Wavelength (nm)", ylab = " ")
#'
#' @export
#'
sim_psr <- function(y, X, x_index = c(1: ncol(X)), nsegs = rep(10, 2), bdegs = rep(3, 3),
                    lambdas = rep(1, 2), pords = rep(2, 2), max_iter = 100)
  {
  # 1: Initialization: estimate alpha (from psSignal) and f (from psNormal)
  psr_fit <- psSignal(y, X, x_index, nseg = nsegs[1], bdeg = bdegs[1],
                      lambda = 1, pord = pords[1], ridge_adj = 0, int = TRUE)
  B <- psr_fit$B
  Q <- X %*% B
  D <- diag(ncol(B))
  D <- diff(D, diff = pords[1])
  alpha <- as.vector(psr_fit$coef[-1])
  yint <- as.vector(psr_fit$coef[1])
  eta <- as.vector(Q %*% alpha) + yint

  f_fit <- psNormal(
    x = eta, y = y, nseg = nsegs[2], bdeg = bdegs[2], pord = pords[2],
    lambda = lambdas[2])

  f_eta = predict(f_fit, x = eta)
  der = psNormal_Deriv(x = eta, y = y, nseg = nsegs[2], bdeg = bdegs[2], pord = pords[2],
                 lambda = lambdas[2])$d_fit


  iter <- 1 # Initialize number of iterations
  cv <- f_fit$cv # Track the cv errors
  d_alpha <- NULL # Track convergence of alpha

  # 2: Start iteration
  for (it in 2:max_iter) {

    ## Update alpha and intercept through lsfit
    y_star <- y - f_eta + der * eta
    Q_tilda <- diag(der) %*% cbind(1, Q)
    P <- cbind(0, sqrt(lambdas[1]) * D)
    nix <- rep(0, nrow(D))
    Q_p <- rbind(Q_tilda, P)
    y_p <- c(y_star, nix)

    fit1 <- lsfit(Q_p, y_p, intercept = FALSE)
    yint_new <- as.vector(fit1$coef[1])
    alpha_new <- as.vector(fit1$coef[2:ncol(Q_p)])

    ## Check the convergence of alpha
    tmp1 <- alpha / sqrt(mean(alpha^2))
    tmp2 <- alpha_new / sqrt(mean(alpha_new^2))
    d_alpha <- c(d_alpha, mean((tmp2 - tmp1)^2) / mean(tmp2^2))

    if (d_alpha[(it - 1)] < 10^(-3)) {
      break
    }

    ## Estimate f through psNormal
    alpha <- alpha_new
    yint <- yint_new
    eta <- as.vector(Q %*% alpha) + yint
    f_fit <- psNormal(
      x = eta, y = y, nseg = nsegs[2], bdeg = bdegs[2], pord = pords[2],
      lambda = lambdas[2]
    )
    f_eta <- predict(f_fit, x = eta)
    der = psNormal_Deriv(x = eta, y = y, nseg = nsegs[2], bdeg = bdegs[2], pord = pords[2],
                        lambda = lambdas[2])$d_fit
    }

  out <- list(y = y, cv = cv, iter = (it - 1), alpha = alpha, yint = yint, B = B,
    f_fit= f_fit, f = f_eta, nsegs = nsegs, bdegs = bdegs, lambdas = lambdas,
    pords = pords, delta_alpha = d_alpha, Q = Q, eta = eta, x_index = x_index)

  class(out) <- "simpsr"
  return(out)
  }
