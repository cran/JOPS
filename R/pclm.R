#' Fit a composite link model
#'
#' @description Fit a smooth latent distribution using the penalized composite link model (PCLM).
#'
#' @param y a vector of counts, length \code{m}.
#' @param C a composition matrix, \code{m} by \code{q}.
#' @param B a B-spline basis matrix, \code{q} by \code{n}.
#' @param lambda the penalty parameter.
#' @param pord the  the order of the difference penalty (default = 2).
#' @param itmax  the maximum number of iterations (default = 50).
#' @param show Set to TRUE or FALSE to display iteration history (default = FALSE).
#'
#' @return
#' A list with the following items:
#' \item{alpha}{the estimated B-spline coefficients, length \code{n}.}
#' \item{gamma}{the estimated latent distribution, length \code{q}.}
#' \item{mu}{estimated values of \code{y}, length \code{m}.}
#' \item{dev}{the deviance of the model.}
#' \item{ed}{the effective model dimension.}
#' \item{aic}{Akaike's Information Criterion.}
#'
#' @details
#' The composite link model assumes that \eqn{E(y) = \mu = C\exp(B \alpha)}, where \eqn{\exp(B\alpha)} is
#' a latent discrete distribution, usually on a finer grid than that for \eqn{y}.
#'
#' Note that \code{sum(gamma) == sum(mu)}.
#'
#' @author Paul Eilers and Jutta Gampe
#'
#' @references
#' Eilers, P. H. C. (2007). III-posed problems with counts, the composite link
#' model and penalized likelihood. \emph{Statistical Modelling}, 7(3), 239–254.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' # Left and right boundaries, and counts, of wide intervals of the data
#' cb <- c( 0, 20, 30, 40, 50, 60)
#' ce <- c(20, 30, 40, 50, 60, 70)
#' y <- c(79, 54, 19, 1, 1, 0)
#'
#' # Construct the composition matrix
#' m <- length(y)
#' n <- max(ce)
#' C <- matrix(0, m, n)
#' for (i in 1:m) C[i, cb[i]:ce[i]] <- 1
#'
#' mids = (cb + ce) / 2 - 0.5
#' widths = ce - cb + 1
#' dens = y / widths / sum(y)
#' x = (1:n) - 0.5
#' B = bbase(x)
#' fit = pclm(y, C, B, lambda = 2, pord = 2, show = TRUE)
#' gamma = fit$gamma / sum(fit$gamma)

#' # Plot density estimate and data
#' plot(x, gamma, type = 'l', lwd = 2, xlab = "Lead Concentration", ylab = "Density")
#' rect(cb, 0, ce, dens, density = rep(10, 6), angle = rep(45, 6))

#' @export

pclm <- function(y, C, B, lambda = 1, pord = 2, itmax = 50,  show = FALSE) {
	# Set up penalty
	n <- dim(B)[2]
	D <- diff(diag(n), diff = pord)
	P <- lambda * t(D) %*% D

	# Initialize
	a <- log(rep(sum(y) / n, n))

	# Iterate
	for (it in 1:itmax) {
		eta <- B %*% a
		gamma <- exp(eta)
		mu  <- as.vector(C %*% gamma)
		M <- diag(1 / mu)
		U <- M %*% C %*% diag(c(gamma)) %*% B
		Q <- t(U) %*% diag(mu) %*% U
		z <- t(U) %*% (y - mu) + Q %*% a
		anew <- solve(Q + P, z)
		dif <- max(abs(anew - a))
		if (show)  cat(it, "  ", dif, "\n")
    a = anew
		if (dif < 1e-7) break
	}

	if (it >= itmax) cat("No convergence after", itmax, "iterations\n")

	# Diagnostics
	H <- solve(Q + P) %*% Q
	ok <- y > 0
	dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
	ed  <- sum(diag(H))
	aic <- dev + 2 * ed

	# Prepare output
	fit = list(alpha = a, gamma = gamma, mu = mu, dev = dev, ed = ed, aic = aic)
	return(fit)
}

