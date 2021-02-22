#' Varying-coefficient single-index signal regression using tensor P-splines.
#'
#' @description \code{sim_vcpsr} is a varying-coefficient single-index
#' signal regression approach that allows both the signal coefficients
#' and the unknown link function to vary with
#' an indexing variable \code{t}, e.g. temperature. Two surfaces
#' are estimated (coefficent and link) that can be sliced at arbitary \code{t}.
#' Anisotripic penalization with P-splines is used on both.
#'
#' @import stats
#'
#' @param y a response vector of length \code{m}, usually continuous.
#' @param X the signal regressors with dimension \code{m} by \code{p1}.
#' @param t_var the varying coeffient indexing variable of length \code{m}.
#' @param x_index an index of length \code{p} for columns of signal matrix;
#' default is simple sequence.
#' @param nsegs a vector of length 4 containing
#' the number of evenly spaced segments between min and max, for each
#' the coefficient surface (row and col) and
#' link surface (row and col), resp. (default \code{rep(10, 4)}.
#' @param bdegs a vector of length 4 containing
#' the degree of B-splines, for each
#' the coefficient surface (row and col) and link surface (row and col), resp.
#' (default cubic \code{rep(3, 4)}).
#' @param lambdas a vector of length 4 containing
#' the positive tuning parameters, for each
#' the coefficient surface (row and col) and link surface (row and col), resp.
#' (default \code{rep(1, 4)}).
#' @param pords a vector of length 4 containing
#' the difference penalty order, for each
#' the coefficient surface (row and col) and link surface (row and col), resp.
#' (default \code{rep(2, 4)}).
#' @param max_iter a scalar for the maximum number of iterations (default 100)
#' @param mins A vector length 2, containing min for signal index and \code{t_var}, default
#' associated with \code{x_index} and \code{t_var} minimums; default is respective minimums.
#' @param maxs A vector length 2, containing max for signal index and \code{t_var}, default
#' associated with \code{x_index} and \code{t_var} maximums; default is respective maximums.

#' @return
#' \item{y}{the response vector of length \code{m}.}
#' \item{alpha}{the P-spline coefficient vector (unfolded) of length \code{(nsegs[1]+bdeg[1])*(negs[2]+bdeg[2])}.}
#' \item{iter}{the number of iterations used for the single-index fit.}
#' \item{yint}{the estimated y-intercept for the single-index model.}
#' \item{Bx}{the B-spline matrix built along the signal index, using \code{nsegs[1]}, used for the coefficient surface.}
#' \item{By}{the B-spline matrix built along the \code{t_var} index,
#' using \code{nsegs[2]}, used for the coefficient surface.}
#' \item{Q}{the effective regressors from the \code{psVCSignal} portion of the single-index
#' fit with dimension \code{m} by \code{length(alpha)}.}
#' \item{t_var}{the VC indexing variable of length \code{m}.}
#' \item{nsegs}{a vector of length 4 containing
#' the number of evenly spaced segments between min and max, for each
#' the coefficient surface (row and col) and
#' link surface (row and col).}
#' \item{bdegs}{a vector of length 4 containing
#' the degree of B-splines, for each
#' the coefficient surface (row and col) and link surface (row and col).}
#' \item{lambdas}{a vector of length 4 containing
#' the positive tuning parameters, for each
#' the coefficient surface (row and col) and link surface (row and col).}
#' \item{pords}{a vector of length 4 containing
#' the difference penalty order, for each
#' the coefficient surface (row and col) and link surface (row and col).}
#' \item{mins}{a vector length 2, containing min for signal index and \code{t_var}.}
#' \item{maxs}{a vector length 2, containing max for signal index and \code{t_var}.}
#' \item{eta}{the estimated linear predictor for the single-index fit.}
#' \item{Pars}{a matrix of 2 rows associated with the signal coefficient surface
#' design parameters, each row: \code{c(min, max, nseg, bdeg, lambda, pord)} for
#' linear predictor \code{x_index} and \code{t_var}, resp.}
#' \item{pPars}{a matrix of 2 rows associated with the link function
#' design parameters, each row: \code{c(min, max, nseg, bdeg, lambda, pord)} for linear
#' predictor \code{eta} and \code{t_var}, resp.}
#' \item{cv}{the leave-one-out cross-validation statistic
#' or the standard error of prediction for the single-index fit.}
#' \item{delta_alpha}{change measure in signal-coefficent parameters at
#' convergence.}
#' \item{fit2D}{\code{ps2DNormal} object, fitting f(\code{eta}, \code{t_var}).}

#' @author Paul Eilers and Brian Marx
#' @references Marx, B. D. (2015). Varying-coefficient single-index signal
#' regression. \emph{Chemometrics and Intelligent Laboratory Systems}, 143, 111–121.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @examples
#' # Load libraries
#' library(fields) # Needed for plotting
#'
#' # Get the data
#' Dat <- Mixture
#'
#' # Dimensions: observations, temperature index, signal
#' m <- 34
#' p1 <- 401
#' p2 <- 12
#'
#' # Stacking mixture data, each mixture has 12 signals stacked
#' # The first differenced spectra are also computed.
#' mixture_data <- matrix(0, nrow = p2 * m, ncol = p1)
#' for (ii in 1:m)
#' {
#'   mixture_data[((ii - 1) * p2 + 1):(ii * p2), 1:p1] <-
#'     t(as.matrix(Dat$xspectra[ii, , ]))
#'   d_mixture_data <- t(diff(t(mixture_data)))
#' }
#'
#' # Response (typo fixed) and index for signal
#' y_mixture <- Dat$fractions
#' y_mixture[17, 3] <- 0.1501
#' index_mixture <- Dat$wl
#'
#' # Select response and replicated for the 12 temps
#' # Column 1: water; 2: ethanediol; 3: amino-1-propanol
#' y <- as.vector(y_mixture[, 2])
#' y <- rep(y, each = p2)
#'
#' bdegs = c(3, 3, 3, 3)
#' pords <- c(2, 2, 2, 2)
#' nsegs <- c(12, 5, 5, 5) # Set to c(27, 7, 7 ,7) for given lambdas
#' mins <- c(700, 30)
#' maxs <- c(1100, 70)
#' lambdas <- c(1e-11, 100, 0.5, 1) # based on svcm search
#' x_index <- seq(from = 701, to = 1100, by = 1) # for dX
#' t_var_sub <- c(30, 35, 37.5, 40, 45, 47.5, 50, 55, 60, 62.5, 65, 70)
#' t_var <- rep(t_var_sub, m)
#' max_iter <- 2 # Set higher in practice, e.g. 100
#' int <- TRUE
#'
#' # Defining x as first differenced spectra, number of channels.
#' x <- d_mixture_data
#'
#'
#' # Single-index VC model using optimal tuning
#' fit <- sim_vcpsr(y, x, t_var, x_index, nsegs, bdegs, lambdas, pords,
#'              max_iter = max_iter, mins = mins, maxs = maxs)
#'
#' plot(fit, xlab = "Wavelength (nm)", ylab = "Temp C")
#' @export
#'
sim_vcpsr <- function(y, X, t_var, x_index = c(1:ncol(X)),
                      nsegs = rep(10, 4), bdegs = rep(3, 4), lambdas = rep(1, 4),
                      pords = rep(2, 4), max_iter = 100,  mins = c(min(x_index), min(t_var)),
                      maxs = c(max(x_index), max(t_var))) {

# 1: Initialization: estimate alpha (from psVCSignal) and
# f (from psNormal)
  psr_fit <- psVCSignal(y, X, x_index, t_var, Pars = rbind(
                  c(mins[1], maxs[1], nsegs[1], bdegs[1], 1, pords[1]),
                  c(mins[2], maxs[2], nsegs[2], bdegs[2], 1, pords[2])),
                  int = TRUE)
  Pars = psr_fit$Pars

  # Effective regressors Q = U%*%B2, with U=X%*%B1
  Q <- psr_fit$Q[, -1]

# Basis along signal and t_var
  Bx <- psr_fit$Bx
  By <- psr_fit$By

# Construction of row and col penalty matrices
  n1 <- nsegs[1] + bdegs[1]
  n2 <- nsegs[2] + bdegs[2]
  d1 <- pords[1]
  D1 <- diag(n1)
  D1 <- diff(D1, diff = d1)
  lambda1 <- lambdas[1]
  P1 <- sqrt(lambda1) * kronecker(D1, diag(n2))
  d2 <- pords[2]
  D2 <- diag(n2)
  D2 <- diff(D2, diff = d2)
  lambda2 <- lambdas[2]
  P2 <- sqrt(lambda2) * kronecker(diag(n1), D2)
  Pen <- rbind(P1, P2)

# Zeroes needed for data augmentation
  z1 <- rep(0, n2 * (n1 - d1))
  z2 <- rep(0, n1 * (n2 - d2))
  nix <- c(z1, z2)

# Coefficients associated with tensor product and y-intercept
  alpha <- as.vector(psr_fit$pcoef[-1])
  yint <- as.vector(psr_fit$pcoef[1])

# Ingredients for pseudo response and scaled basis
  eta <- as.vector(Q %*% alpha) + yint
  xxyy <- cbind(eta, t_var)
  mineta <- min(eta)
  maxeta <- max(eta)
  pPars <- rbind(c(mineta, maxeta, nsegs[3], bdegs[3], lambdas[3], pords[3]),
                 c(min(t_var), max(t_var), nsegs[4], bdegs[4], lambdas[4], pords[4]))
  teta <- ps2D_PartialDeriv(Data = cbind(xxyy, y), Pars = pPars, xxyy)
  f_fit <- teta
  der <- as.vector(teta$d_fit)
  f_eta <- teta$fit

# Initialize number of iterations
 iter <- 1
#Track the cv errors
  cv <- f_fit$cv
# Track convergence of alpha
  d_alpha <- NULL

# 2: Start iteration
  for (it in 2:max_iter) {
# Update alpha and intercept through lsfit
    y_star <- y - f_eta + der * eta
    Q_tilda <- diag(der) %*% cbind(1, Q)
    p <- cbind(0, Pen)
    Q_p <- rbind(Q_tilda, p)
    y_p <- c(y_star, nix)
    fit1 <- lsfit(Q_p, y_p, intercept = FALSE)
    yint_new <- as.vector(fit1$coef[1])
    alpha_new <- as.vector(fit1$coef[2:ncol(Q_p)])

# Check the convergence of alpha
    tmp1 <- alpha / sqrt(mean(alpha^2))
    tmp2 <- alpha_new / sqrt(mean(alpha_new^2))
    d_alpha <- c(d_alpha, mean((tmp2 - tmp1)^2) / mean(tmp2^2))
    if (d_alpha[(it - 1)] < 10^(-3)) {
      break
    }

# Estimate f through pnormal
    alpha <- alpha_new
    yint <- yint_new
    eta <- as.vector(Q %*% alpha) + yint
    mineta <- min(eta)
    maxeta <- max(eta)
    mintvar <- min(t_var)
    maxtvar <- max(t_var)
    xxyy <- cbind(eta, t_var)
    pPars <- rbind(
      c(mineta, maxeta, nsegs[3], bdegs[3], lambdas[3], pords[3]),
      c(mintvar, maxtvar, nsegs[4], bdegs[4], lambdas[4], pords[4])
    )
    teta <- ps2D_PartialDeriv(
              Data = cbind(xxyy, y), Pars = pPars, xxyy)
    f_fit <- teta
    fit2D = ps2DNormal(Data = cbind(xxyy, y), Pars = pPars)
    der <- as.vector(teta$d_fit)
    f_eta <- teta$fit
  }
  out <- list(
    y = y, iter = (it - 1), alpha = alpha, yint = yint, Bx = Bx,
    By = By, Q = Q,  nsegs = nsegs, lambdas = lambdas, pords = pords,
    t_var = t_var, eta = eta, Pars = Pars, pPars = pPars, bdegs = bdegs,
    mins = mins, maxs = maxs, cv = cv, delta_alpha = d_alpha, fit2D = fit2D)

  class(out) <- "simvcpsr"
  return(out)
}


