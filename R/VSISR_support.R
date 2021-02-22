#' Partial derivative two-dimensional smoothing scattered (normal)
#' data using P-splines.
#'
#' @description \code{ps2D_PartialDeriv} provides the partial derivative
#' P-spline surface along \code{x}, with aniosotripic penalization of
#' tensor product B-splines.
#'
#' @details This is support function for \code{sim_vcpsr}.
#'
#' @param Data a matrix of 3 columns \code{x, y, z} of equal length;
#' the response is \code{z}.
#' @param Pars a matrix of 2 rows, where the first and second row
#' sets the P-spline paramters for \code{x} and \code{y}, respectively.
#' Each row consists of: \code{min max nseg bdeg lambda pord}.
#' The \code{min} and \code{max} set the ranges, \code{nseg} (default 10)
#' is the number of evenly spaced segments between \code{min} and \code{max},
#' \code{bdeg} is the degree of the basis (default 3 for cubic),
#' \code{lambda} is the (positive) tuning parameter for the penalty (default 1),
#' \code{pord} is the number for the order of the difference penalty (default 2).
#' @param XYpred a matrix with two columns \code{(x, y)} that give the coordinates
#' of (future) prediction; the default is the data locations.

#' @return
#' \item{coef}{a vector of length \code{(Pars[1, 3] + Pars[1, 4]) * (Pars[1, 3] + Pars[1, 4]).}
#' of (unfolded) estimated P-spline coefficients.}
#' \item{B}{the tensor product B-spline matrix of dimensions \code{m} by \code{length(coef)}.}
#' \item{fit}{a vector of \code{length(y)} of smooth estimated means (at the \code{x, y} locations).}
#' \item{pred}{a vector of length \code{nrow(XYpred)} of (future) predictions.}
#' \item{d_coef}{a vector of length \code{(Pars[1, 3] + Pars[1,4] - 1) * (Pars[1,3]+Pars[1,4]).}
#' of (unfolded) partial derivative estimated P-spline coefficients.}
#' \item{B_d}{the tensor product B-spline matrix of dimensions \code{m} by \code{lengh(d_coef)}, associated with
#' the partial derivative of the tensor basis.}
#' \item{d_fit}{a vector of \code{length(y)} of partial derivative (along \code{x})
#' of the smooth estimated means (at the \code{x, y} locations).}
#' \item{d_pred}{a vector of length \code{nrow(XYpred)} of partial derivative (future) predictions.}
#' \item{Pars}{a matrix of 2 rows, where each the first (second) row
#'  sets the P-spline paramters for \code{x (y)}: \code{min max nseg bdeg lambda pord}. See the argument above.}
#' \item{cv}{root leave-one-out CV or root average PRESS.}
#' \item{XYpred}{a matrix with two columns \code{(x, y)} that give the coordinates
#' of (future) prediction; the default is the data locations.}

#' @author Brian Marx

#' @references Marx, B. D. (2015). Varying-coefficient single-index signal
#' regression. \emph{Chemometrics and Intelligent Laboratory Systems}, 143, 111–121.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.

#' @export
"ps2D_PartialDeriv" <- function(Data, Pars = rbind(c(min(Data[, 1]), max(Data[, 1]), 10, 3, 1, 2),
                                                   c(min(Data[, 2]), max(Data[, 2]), 10, 3, 1, 2)),
                                                   XYpred = cbind(Data[, 1], Data[, 2])) {
  # P-spline tensor product fit on (x, y, z)
  fit2D <- ps2DNormal(Data = Data, Pars = Pars)
  pcoef <- fit2D$coef
  pfit <- fit2D$fit
  cv <- fit2D$cv
  zpred <- fit2D$pred
  B <- fit2D$B
  nx <- fit2D$Pars[1, 3] + fit2D$Pars[1, 4]
  ny <- fit2D$Pars[2, 3] + fit2D$Pars[2, 4]

  # Prepare bases and get fitted values for XYpred surfaces
  Bxp <- bbase(XYpred[, 1], Pars[1, 1], Pars[1, 2], Pars[1, 3], Pars[1, 4])
  Byp <- bbase(XYpred[, 2], Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
  B1p <- kronecker(Bxp, t(rep(1, ny)))
  B2p <- kronecker(t(rep(1, nx)), Byp)
  Bp <- B1p * B2p
  zpred <- Bp %*% pcoef

  # Tensor product basis for fitted partial derivative (along x) surface.
  Pars1_d <- Pars[1, ]
  Pars1_d[4] <- Pars1_d[4] - 1
  Bx_d <- bbase(Data[, 1], Pars1_d[1], Pars1_d[2], Pars1_d[3], Pars1_d[4])
  By <- bbase(Data[, 2], Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
  nx_d <- ncol(Bx_d)
  ny <- ncol(By)
  B1d <- kronecker(Bx_d, t(rep(1, ny)))
  B2d <- kronecker(t(rep(1, nx_d)), By)
  B_d <- B1d * B2d

  # Array of coefficents for partial derivative surface
  pcoef_m <- matrix(pcoef, nx, ny, byrow = TRUE)
  d_pcoef_m <- (diff((pcoef_m)))
  d_pcoef <- c(t(d_pcoef_m)) / ((Pars[1, 2] - Pars[1, 1]) / (1 * Pars[1, 3]))
  d_fit <- as.vector(B_d %*% d_pcoef)

  # Tensor product basis for XYpred partial derivative (along x) surface.
  Bxp_d <- bbase(XYpred[, 1], Pars1_d[1], Pars1_d[2], Pars1_d[3], Pars1_d[4])
  Byp <- bbase(XYpred[, 2], Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
  B1p_d <- kronecker(Bxp_d, t(rep(1, ny)))
  B2p <- kronecker(t(rep(1, nx_d)), Byp)
  Bp_d <- B1p_d * B2p
  zpred_d <- as.vector(Bp_d %*% d_pcoef)
  P <- list(
    coef = pcoef, B = B, fit = pfit, d_coef = d_pcoef, B_d = B_d,
    d_fit = d_fit, pred = zpred, d_pred = zpred_d, cv = cv, XYpred = XYpred,
    Pars = Pars
  )
  P
}


#' Derivative for a P-spline fit of scattered (normal)
#' data.
#' @description \code{psNormal_Deriv} provides the derivative
#' P-spline fit along \code{x}.
#'
#' @details This is also a
#' support function needed for \code{sim_psr} and \code{sim_vcpsr}.
#' SISR (Eilers, Li, Marx, 2009).
#'
#' @param y the response vector, usually continuous data.
#' @param x the vector for the continuous regressor of \code{length(y)} and the abcissae of fit.
#' @param wts the vector of weights, default is 1; 0/1 allowed.
#' @param xl the number for the min along \code{x} (default is min(\code{x})) .
#' @param xr the number for the max along \code{x} (default is max(\code{x})).
#' @param nseg the number of evenly spaced segments between \code{xl} and \code{xr}.
#' @param bdeg the number of the degree of the basis, usually 1, 2, or 3 (defalult).
#' @param pord the number of the order of the difference penalty, usually 1, 2 (defalult), or 3.
#' @param lambda the positive tuning parameter (default 1).
#' @param xgrid a scalar or a vector that gives the \code{x} locations for prediction, useful for plotting.
#' If a scalar (default 100) is used then a uniform grid of this size along (\code{xl}, \code{xr}).

#' @return
#' \item{coef}{a vector of \code{length(nsegs + bdeg)}
#' of estimated P-spline coefficients.}
#' \item{B}{The B-spline matrix of dimensions \code{m} by \code{length(coef)}.}
#' \item{fit}{a vector of \code{length(y)} of smooth estimated means (at the \code{x} locations).}
#' \item{pred}{a vector of \code{length(xgrid)} of (future) predictions.}
#' \item{d_coef}{a vector of \code{length(nsegs + bdeg - 1)}
#' of differenced (derivative) estimated P-spline coefficients.}
#' \item{B_d}{The first derivative B-spline matrix of dimensions \code{m} by \code{lengh(d_coef)}.}
#' \item{d_fit}{a vector of \code{length(y)} of partial derivative (along \code{x})
#' of the smooth estimated means (at the \code{x} locations).}
#' \item{d_pred}{a vector of length \code{lenght(xgrid)} of partial derivative (future) predictions.}
#' \item{xl}{the number for the min along \code{x} (default is min(\code{x})).}
#' \item{xr}{the number for the max along \code{x} (default is max(\code{x})).}
#' \item{nseg}{the number of evenly spaced segments between \code{xl} and \code{xr}.}
#' \item{bdeg}{the number of the degree of the basis, usually 1, 2, or 3 (default).}
#' \item{pord}{the number of the order of the difference penalty, usually 1, 2 (default), or 3.}
#' \item{lambda}{the positive tuning parameter (default 1).}


#' @author Paul Eilers and Brian Marx
#' @seealso sim_psr sim_vcpsr
#' @references Marx, B. D. (2015). Varying-coefficient single-index signal
#' regression. \emph{Chemometrics and Intelligent Laboratory Systems}, 143, 111–121.
#' @references Eilers, P.H.C., B. Li, B.D. Marx (2009).
#' Multivariate calibration with single-index signal regression,
#' \emph{Chemometrics and Intellegent Laboratory Systems}, 96(2), 196-202.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.

#' @export
"psNormal_Deriv" <-
  function(x, y, xl = min(x), xr = max(x), nseg = 10, bdeg = 3, pord = 2,
             lambda = 1, wts = rep(1, length(y)), xgrid = x) {
    f <- psNormal(x = x, y = y, nseg = nseg, bdeg = bdeg, pord = pord, lambda = lambda, wts = wts)
    m <- length(y)
    alpha <- as.vector(f$pcoef)
    mu <- f$muhat
    cv <- f$cv
    fit_grid <- f$ygrid

    # First derivative prediction at x
    B_der <- bbase(x, xl, xr, nseg, bdeg - 1)
    alpha_der <- diff(alpha) / ((xr - xl) / nseg)
    d_fit <- as.vector(B_der %*% alpha_der)

    B_pder <- bbase(xgrid, xl, xr, nseg, bdeg - 1)
    d_pfit <- as.vector(B_pder %*% alpha_der)

    # Return list
    P <- list(
      coef = alpha, B = f$B, fit = mu, d_coef = alpha_der,
      B_d = B_der, xl = xl, xr = xr,
      d_fit = d_fit, pred = fit_grid,
      d_pred = d_pfit, cv = cv, xgrid = xgrid,
      nseg = nseg, bdeg = bdeg, pord = pord, lambda = lambda
    )
    P
  }
