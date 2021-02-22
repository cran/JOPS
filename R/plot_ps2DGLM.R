#' Plotting function for \code{ps2DGLM}
#'
#' @description Plotting function for 2D P-spline (GLM) smooothing
#' (using \code{ps2DGLM} with \code{class ps2dglm}).
#'
#' @import graphics
#'
#' @param x the P-spline object, usually from \code{ps2DGLM}.
#' @param ... other parameters.
#' @param xlab label for the x-axis, e.g. "my x" (quotes required).
#' @param ylab label for the y-axis, e.g. "my y" (quotes required).
#' @param se a scalar, e.g. \code{se = 2} to produce twice se surfaces,
#' set \code{se} > 0 (or set \code{se = 0} to supress).
#' @param Resol resolution for plotting, default \code{Resol = 100}.

#' @return
#' \item{Plot}{a plot of the mean (inverse link) 2D P-spline (GLM) smooth surface.}
#'
#' @author Paul Eilers and Brian Marx
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#'
#' @examples
#' library(fields)
#' library(JOPS)
#' # Extract data
#' library(rpart)
#' Kyphosis <- kyphosis$Kyphosis
#' Age <- kyphosis$Age
#' Start <- kyphosis$Start
#' y <- 1 * (Kyphosis == "present") # make y 0/1
#' fit <- ps2DGLM(
#'   Data = cbind(Start, Age, y),
#'   Pars = rbind(c(1, 18, 10, 3, .1, 2), c(1, 206, 10, 3, .1, 2)),
#'   family = "binomial"
#' )
#' plot(fit, xlab = "Start", ylab = "Age")
#' #title(main = "Probability of Kyphosis")
#' @export

# library(fields)
plot.ps2dglm <- function(x,..., xlab = " ", ylab = " ", Resol = 100, se = 2) {
  ps2dg = x
  Pars <- ps2dg$Pars
  Data <- ps2dg$Data
  x <- Data[, 1]
  y <- Data[, 2]
  z <- Data[, 3]
  qr <- ps2dg$qr
  Q = ps2dg$Q
  dispersion_parm <- ps2dg$dispersion_parm
  n_col <- ncol(Q)
  pcoef <- ps2dg$pcoef
  link <- ps2dg$link

  # Compute grid, useful for plotting
  u <- seq(min(x), max(x), length = Resol)
  v <- seq(min(y), max(y), length = Resol)
  U_ <- outer(rep(1, Resol), u)
  V_ <- outer(v, rep(1, Resol))
  U <- as.vector(U_)
  V <- as.vector(V_)
  Bxgrid <- bbase(U, Pars[1, 1], Pars[1, 2], Pars[1, 3], Pars[1, 4])
  Bygrid <- bbase(V, Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
  n1 <- ncol(Bxgrid)
  n2 <- ncol(Bygrid)
  B1grid <- kronecker(Bxgrid, t(rep(1, n2)))
  B2grid <- kronecker(t(rep(1, n1)), Bygrid)
  Bgrid <- B1grid * B2grid
  zgrid <- as.vector(Bgrid %*% pcoef)
  izgrid <- inverse_link(x = zgrid, link = link)
  Fitgrid <- matrix(izgrid, Resol, Resol, byrow = TRUE)

  # Plot image of inverse link (predicted)
  image.plot(U_[1, ], V_[, 1], Fitgrid,
    xlab = xlab, ylab = ylab,
    sub = "Fitted"
  )

  # Lower and Upper se*SE surfaces
  if (se > 0) {
    R <- qr.R(qr)
    # C2 = chol2inv(R)

    # Variances of fitted values
    L <- forwardsolve(t(R), t(Bgrid))
    var_hat <- colSums(L * L)
    stdev_hat <- sqrt(dispersion_parm) * as.vector(sqrt(var_hat))
    pivot <- se * stdev_hat
    upper <- zgrid + pivot
    lower <- zgrid - pivot
    ilower <- inverse_link(x = lower, link = link)
    iupper <- inverse_link(x = upper, link = link)

    if (link == "reciprocal") {
      ilowup <- cbind(ilower, iupper)
      iupper <- ilowup[, 1]
      ilower <- ilowup[, 2]
    }
    L_hatm <- matrix(ilower, Resol, Resol, byrow = TRUE)
    U_hatm <- matrix(iupper, Resol, Resol, byrow = TRUE)

    image.plot(U_[1, ], V_[, 1], U_hatm,
      xlab = xlab,
      ylab = ylab, sub = "2 se Upper Surface"
    )

    image.plot(U_[1, ], V_[, 1], L_hatm,
      xlab = xlab,
      ylab = ylab, sub = "2 se Lower Surface"
    )
  }
  if (se < 0) {
    warning(paste("se should be nonnegative"))
  }
}
