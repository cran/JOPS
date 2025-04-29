#' Two-dimensional smoothing scattered (normal) data using P-splines.
#'
#' @description ps2DNormal is used to smooth scattered
#' (normal) data, with anisotropic penalization of
#' tensor product P-splines.
#'
#'
#' @details Support functions needed: \code{pspline_fitter}, \code{bbase}, and \code{pspline_2dchecker}.
#' @seealso ps2DGLM
#'
#' @import stats
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
#' \code{pord} is the number for the order of the difference penalty (default 2),
#' @param XYpred a matrix with two columns \code{(x,y)} that give the coordinates
#' of (future) prediction; the default is the data locations.

#' @return
#' \item{coef}{a vector of length \code{(Pars[1,3]+Pars[1,4])*(Pars[2,3]+Pars[2,4])}
#' of (unfolded) estimated P-spline coefficients.}
#' \item{fit}{a vector of \code{length(y)} of smooth estimated means (at the \code{x,y} locations).}
#' \item{pred}{a vector of length \code{nrow(XYpred)} of (future) predictions.}
#' \item{Pars}{the design and tuning parameters (see arguments above).}
#' \item{cv}{leave-one-out standard error of prediction or root average PRESS.}
#' \item{h}{"hat" diagonals of tensor P-spline fit.}
#' \item{B}{tensor product B-spline basis used for fitting.}
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#' @import stats
#'
#' @examples
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
#'
ps2DNormal <- function(Data, Pars = rbind(c(min(Data[, 1]), max(Data[, 1]), 10, 3, 1, 2),
                                          c(min(Data[, 2]), max(Data[, 2]), 10, 3, 1, 2)),
                                          XYpred=expand.grid(Data[, 1],Data[, 2])) {
  # Prepare bases
  p1 <- Pars[1, ]
  p2 <- Pars[2, ]
  Bx <- bbase(Data[, 1], p1[1], p1[2], p1[3], p1[4])
  By <- bbase(Data[, 2], p2[1], p2[2], p2[3], p2[4])
  m = nrow(Bx)
  nx <- ncol(Bx)
  ny <- ncol(By)

  # Compute tensor products
  B1 <- kronecker(Bx, t(rep(1, ny)))
  B2 <- kronecker(t(rep(1, nx)), By)
  B <- B1 * B2

  # Construct penalty matrices
  dx <- Pars[1, 6]
  Dx <- diag(nx)
  Dx <- diff(Dx, diff = dx)
  lambdax <- Pars[1, 5]
  Px <- sqrt(lambdax) * kronecker(Dx, diag(ny))
  dy <- Pars[2, 6]
  Dy <- diag(ny)
  Dy <- diff(Dy, diff = dy)

  lambday <- Pars[2, 5]
  Py <- sqrt(lambday) * kronecker(diag(nx), Dy)

  # Data augmentation and regression
  zx <- rep(0, ny * (nx - dx))
  zy <- rep(0, nx * (ny - dy))
  zplus <- c(Data[, 3], zx, zy)
  Bplus <- rbind(B, Px, Py)
  fit1 <- lsfit(Bplus, zplus, intercept = FALSE)
  pcoef = fit1$coef
  pfit <- B%*%pcoef
  h <- hat(fit1$qr, intercept = FALSE)[1:m]
  press_e <- fit1$residuals[1:m]/(1 - h)
  cv <- sqrt(sum((press_e)^2)/(m))

  # Prediction
  Bxp <- bbase(XYpred[, 1], p1[1], p1[2], p1[3], p1[4])
  Byp <- bbase(XYpred[, 2], p2[1], p2[2], p2[3], p2[4])
  B1p <- kronecker(Bxp, t(rep(1, ny)))
  B2p <- kronecker(t(rep(1, nx)), Byp)
  Bp <- B1p * B2p
  zpred <- Bp %*% pcoef

  P <- list(
    coef = pcoef, fit = pfit, pred = zpred,
    XYpred = XYpred, Pars = Pars, Data = Data, cv = cv,
    h = h, B = B
  )
  class(P) <- "ps2dnormal"
  return(P)
}
