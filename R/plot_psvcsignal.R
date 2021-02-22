#' Plotting function for \code{psVCSignal}
#'
#' @description Plotting function for varying-coefficent signal
#' regression P-spline smooth coefficients (using \code{psVCSignal} with \code{class psvcsignal}).
#' Although se surface bands can be comuputed they are intentially left out as they are not
#' interpretable, and there is generally little data to steer
#' such a high-dimensional parameterization.
#'
#' @import graphics
#'
#' @param x the P-spline object, usually from \code{psVCSignal}.
#' @param ... other parameters.
#' @param xlab label for the x-axis, e.g. "my x" (quotes required).
#' @param ylab label for the y-axis, e.g. "my y" (quotes required).
#' @param Resol resolution for plotting, default \code{Resol = 100}.


#' @return
#' \item{Plot}{a two panel plot, one of the 2D P-spline signal coefficient surface
#' and another that displays several slices of the smooth coefficient vectors at fixed levels of the
#' varying index.}
#'
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P. H. C. and Marx, B. D. (2003). Multivariate calibration with temperature
#' interaction using two-dimensional penalized signal regression. \emph{Chemometrics
#' and Intellegent Laboratory Systems}, 66, 159–174.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' library(fds)
#' data(nirc)
#' iindex <- nirc$x
#' X <- nirc$y
#' sel <- 50:650 # 1200 <= x & x<= 2400
#' X <- X[sel, ]
#' iindex <- iindex[sel]
#' dX <- diff(X)
#' diindex <- iindex[-1]
#' y <- as.vector(labc[1, 1:40]) # percent fat
#' t_var <- as.vector(labc[4, 1:40]) # percent flour
#' oout <- 23
#' dX <- t(dX[, -oout])
#' y <- y[-oout]
#' t_var = t_var[-oout]
#' Pars = rbind(c(min(diindex), max(diindex), 25, 3, 1e-7, 2),
#' c(min(t_var), max(t_var), 20, 3, 0.0001, 2))
#' fit1 <- psVCSignal(y, dX, diindex, t_var, Pars = Pars,
#' family = "gaussian", link = "identity", int = TRUE)
#' plot(fit1, xlab = "Coefficient Index", ylab = "VC: % Flour")
#' names(fit1)
#' @export
#'

#library(fields)
plot.psvcsignal = function(x, ..., xlab = " ", ylab = " ",
                           Resol = 100){

# Prepare bases for estimated alpha surface
psvcsig = x
Pars = psvcsig$Pars
pcoef = psvcsig$pcoef
int = psvcsig$int
Q = psvcsig$Q
n_col=ncol(Q)

S_index_ <- seq(from = Pars[1, 1], to = Pars[1, 2], length = Resol)
oS <- outer(rep(1, Resol), S_index_)
Bx_ <- bbase(as.vector(oS), Pars[1, 1],Pars[1, 2],Pars[1, 3], Pars[1, 4])
t_index_ <- seq(from = Pars[2, 1], to = Pars[2, 2], length = Resol)
ot <- outer(t_index_, rep(1, Resol))
By_ <- bbase(as.vector(ot), Pars[2, 1 ],Pars[2, 2],Pars[2, 3], Pars[2, 4])

# Compute tensor products for estimated alpha surface
n1=ncol(Bx_)
n2=ncol(By_)
B1_ <- kronecker(Bx_, t(rep(1, n2)))
B2_ <- kronecker(t(rep(1, n1)), By_)
B_ <- B1_ * B2_
    if (int) {
      A_hat <- as.vector(B_ %*% pcoef[(1+int):n_col])
    }
    if (!int) {
      A_hat <- as.vector(B_ %*% pcoef)
    }

oldpar = par(no.readonly = TRUE)
on.exit(par(oldpar))
par(mfrow = c(1,2))
# Put coefficients in matrix and create image plot
A_hatm <- matrix(A_hat, Resol, Resol, byrow = TRUE)
image.plot(oS[1,  ], ot[, 1], A_hatm, xlab = xlab, ylab = ylab)

matplot(S_index_, A_hatm[, seq(1, Resol, length = 6)],
        type = "l", col = (1:6), lty = c(1:6), ylab = " ", xlab = xlab
)
abline(0, 0)

}
