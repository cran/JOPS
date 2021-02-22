#' Plotting function for \code{sim_vcpsr}
#'
#' @description Plotting function for varying-coefficent single-index signal
#' regression using tensor P-splines (using \code{sim_vcpsr} with \code{class simvcpsr}).
#'
#' @import graphics
#'
#' @param x the P-spline object, usually from \code{sim_vcpsr}.
#' @param ... other parameters.
#' @param xlab label for the x-axis, e.g. "my x" (quotes required).
#' @param ylab label for the y-axis, e.g. "my y" (quotes required).
#' @param Resol resolution for plotting, default \code{Resol = 100}.
#'
#' @return
#' \item{Plot}{a plot of the estimated 2D P-spline signal coefficent surface along with the companion plot of the estimated
#' 2D P-spline varying link function surface. Slices of these plots, at fixed levels of the indexing covariate, are also provided.}
#'
#' @author Paul Eilers and Brian Marx
#' @references Marx, B. D. (2015). Varying-coefficient single-index signal
#' regression. \emph{Chemometrics and Intellegent Laboratory Systems}, 143, 111–121.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' #' @examples
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

#'
#'
#' @export
#'
# library(fields)
plot.simvcpsr <- function(x,..., xlab = " ", ylab = " ",
                          Resol = 100) {
  simvc = x
  eta <- simvc$eta
  nsegs <- simvc$nsegs
  mins <- simvc$mins
  maxs <- simvc$maxs
  pords <- simvc$pords
  lambdas <- simvc$lambdas
  y <- simvc$y
  bdegs <- simvc$bdegs
  alpha <- simvc$alpha
  t_var = simvc$t_var

  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(2, 2))
  Pars1 <- c(mins[1], maxs[1], nsegs[1], bdegs[1])
  Pars2 <- c(mins[2], maxs[2], nsegs[2], bdegs[2])

  # Prepare bases for estimated alpha surface
  S_index_ <- seq(from = Pars1[1], to = Pars1[2], length = Resol)
  oS <- outer(rep(1, Resol), S_index_)
  Bx_ <- bbase(as.vector(oS), Pars1[1], Pars1[2], Pars1[3], Pars1[4])
  t_index_ <- seq(from = Pars2[1], to = Pars2[2], length = Resol)
  ot <- outer(t_index_, rep(1, Resol))
  By_ <- bbase(as.vector(ot), Pars2[1], Pars2[2], Pars2[3], Pars2[4])

  # Compute tensor products for estimated alpha surface
  B1_ <- kronecker(Bx_, t(rep(1, ncol(By_))))
  B2_ <- kronecker(t(rep(1, ncol(Bx_))), By_)
  B_ <- B1_ * B2_
  A_hat <- B_ %*% alpha
  A_hatm <- matrix(A_hat, Resol, Resol, byrow =  TRUE)

  image.plot(oS[1, ], ot[, 1], A_hatm,
    xlab = xlab, ylab = ylab,
    main = "Coefficient surface"
  )

  matplot(seq(Pars1[1], Pars1[2], length = Resol), A_hatm[, seq(1, Resol, length = 6)],
    type = "l", col = (1:6), lty = c(1:6), ylab = " ", xlab = xlab,
    main = "Coef slices, by indexing variable"
  )
  abline(0, 0)

  pPars <- rbind(
    c(min(eta), max(eta), nsegs[3], bdegs[3], lambdas[3], pords[3]),
    c(min(t_var), max(t_var), nsegs[4], bdegs[4], lambdas[4], pords[4])
  )
  eta_index_ <- seq(from = pPars[1, 1], to = pPars[1, 2], length = Resol)
  oeta <- outer(rep(1, Resol), eta_index_)
  t_index_ <- seq(from = pPars[2, 1], to = pPars[2, 2], length = Resol)
  ot <- outer(t_index_, rep(1, Resol))
  teta <- ps2DNormal(cbind(eta, t_var, y),
    Pars = pPars,
    XYpred = cbind(as.vector(oeta), as.vector(ot))
  )
  fit_matrix <- matrix(teta$pred, Resol, Resol, byrow = TRUE)
  image.plot(oeta[1, ], ot[, 1], fit_matrix, xlab = "Linear predictor", ylab = ylab, main = "Link surface")
  matplot(eta_index_, fit_matrix[, seq(1, Resol, length = 6)],
    type = "l", col = (1:6),
    lty = c(1:6), ylab = " ", xlab = "Linear predictor", main = "Link slices"
  )
  abline(0, 1, lty = 2, col = 1, lwd = 2)
}
