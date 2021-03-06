% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pclm.R
\name{pclm}
\alias{pclm}
\title{Fit a composite link model}
\usage{
pclm(y, C, B, lambda = 1, pord = 2, itmax = 50, show = FALSE)
}
\arguments{
\item{y}{a vector of counts, length \code{m}.}

\item{C}{a composition matrix, \code{m} by \code{q}.}

\item{B}{a B-spline basis matrix, \code{q} by \code{n}.}

\item{lambda}{the penalty parameter.}

\item{pord}{the  the order of the difference penalty (default = 2).}

\item{itmax}{the maximum number of iterations (default = 50).}

\item{show}{Set to TRUE or FALSE to display iteration history (default = FALSE).}
}
\value{
A list with the following items:
\item{alpha}{the estimated B-spline coefficients, length \code{n}.}
\item{gamma}{the estimated latent distribution, length \code{q}.}
\item{mu}{estimated values of \code{y}, length \code{m}.}
\item{dev}{the deviance of the model.}
\item{ed}{the effective model dimension.}
\item{aic}{Akaike's Information Criterion.}
}
\description{
Fit a smooth latent distribution using the penalized composite link model (PCLM).
}
\details{
The composite link model assumes that \eqn{E(y) = \mu = C\exp(B \alpha)}, where \eqn{\exp(B\alpha)} is
a latent discrete distribution, usually on a finer grid than that for \eqn{y}.

Note that \code{sum(gamma) == sum(mu)}.
}
\examples{
# Left and right boundaries, and counts, of wide intervals of the data
cb <- c( 0, 20, 30, 40, 50, 60)
ce <- c(20, 30, 40, 50, 60, 70)
y <- c(79, 54, 19, 1, 1, 0)

# Construct the composition matrix
m <- length(y)
n <- max(ce)
C <- matrix(0, m, n)
for (i in 1:m) C[i, cb[i]:ce[i]] <- 1

mids = (cb + ce) / 2 - 0.5
widths = ce - cb + 1
dens = y / widths / sum(y)
x = (1:n) - 0.5
B = bbase(x)
fit = pclm(y, C, B, lambda = 2, pord = 2, show = TRUE)
gamma = fit$gamma / sum(fit$gamma)
# Plot density estimate and data
plot(x, gamma, type = 'l', lwd = 2, xlab = "Lead Concentration", ylab = "Density")
rect(cb, 0, ce, dens, density = rep(10, 6), angle = rep(45, 6))
}
\references{
Eilers, P. H. C. (2007). III-posed problems with counts, the composite link
model and penalized likelihood. \emph{Statistical Modelling}, 7(3), 239–254.

Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
P-splines.} Cambridge University Press.
}
\author{
Paul Eilers and Jutta Gampe
}
