#' Smooth a 2D histogram
#'
#' @description Fit a 2D smooth P-spline surface to a matrix of counts, assuming Poisson distributed observations.
#'
#' @param Y a matrix of counts.
#' @param nsegx the number of knots along \code{x} (default=10).
#' @param nsegy the number of evenly spaced knots along \code{y} for Tensor product B-spline basis (default=10).
#' @param lambdax the positive number for the tuning parameter along \code{x}.
#' @param lambday the positive number for the tuning parameter along \code{y}.
#' @param bdeg the degree of the basis, default is 3.
#' @param dx the order of the difference penalty along \code{x}, default is 3.
#' @param dy the order of the difference penalty along \code{y}, default is 3.
#' @param tol the convergence criterion (default \code{1e-5}).
#' @param Mu the initialization of the mean (default \code{Y + 0.01}).
#' @param kappa a (small, positive) number for ridge tuning parameter to stabilize estimation (default \code{1e-4}).
#'
#' @return A list with elements:
#' \item{ed}{the effective dimension of the smooth 2D surface.}
#' \item{Mu}{a matrix with the smooth estimates, with dimensions of \code{dim(Y)}}
#' \item{pen}{the numerical value of the penalty.}
#'
#' @author Paul Eilers
#'
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' x = faithful$eruptions
#' y = faithful$waiting
#' h = hist2d(x, y, c(100, 100))
#' sm = hist2dsm(h$H, nsegx = 25, nsegy = 25, bdeg = 3, lambdax = 10, lambday = 10)
#' image(h$xgrid, h$ygrid, sm$Mu, xlab = 'Eruption length (min)',
#'       ylab = 'Waiting time (min)', main = 'Old Faithful')
#'
#' @export

hist2dsm = function(Y, nsegx = 10, nsegy = nsegx, bdeg = 3,  lambdax = 10, lambday = lambdax,
                    dx = 3, dy = dx, Mu = Y + 0.01, kappa = 1e-4, tol = 1e-5) {
    nx = nrow(Y)
    ny = ncol(Y)

    # Conmpute the basis matrices
    Bx = bbase(1:nx, 0, nx + 1, nsegx, bdeg)
    By = bbase(1:ny, 0, ny + 1, nsegy, bdeg)
    nbx = ncol(Bx)
    nby = ncol(By)
    Tx = rowtens(Bx)
    Ty = rowtens(By)

    # Compute the penalty matrices
    Dx = diff(diag(nbx), diff = dx)
    Dy = diff(diag(nby), diff = dy)
    Px = lambdax * t(Dx) %*% Dx
    Py = lambday * t(Dy) %*% Dy
    P = kronecker(Py, diag(nbx)) + kronecker(diag(nby), Px)
    P = P + kappa * diag(nrow(P))

    # Initialize
    Z = log(Mu)
    Z = Z - log(sum(exp(Z)) / sum(Y))

    # Iterate, using the array algorithm
    for (it in 1:30) {
      Mu = exp(Z)
      U = Y - Mu + Mu * Z
      Q = t(Tx) %*% Mu %*% Ty
      dim(Q) = c(nbx, nbx, nby, nby)
      Q = aperm(Q, c(1, 3, 2, 4))
      dim(Q) = c(nbx * nby, nbx * nby)
      r = t(Bx) %*% U %*% By
      dim(r) = c(nbx * nby, 1)
      A = solve(Q + P, r)
      a = A
      dim(A) = c(nbx, nby)
      Znew = Bx %*% A %*% t(By)
      dz = sum(abs(Z - Znew))
      if (dz < tol) break
      Z = Znew
    }

    # Compute diagnostics
    K = solve(Q + P, Q)
    ed = sum(diag(K))
    pen = t(a) %*% P %*% a
    return(list(ed = ed, Mu = Mu, pen = pen))
}

