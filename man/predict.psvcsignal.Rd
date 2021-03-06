% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_support.R
\name{predict.psvcsignal}
\alias{predict.psvcsignal}
\title{Predict function for \code{psVCSignal}}
\usage{
\method{predict}{psvcsignal}(object, ..., X_pred, t_pred, type = "mu")
}
\arguments{
\item{object}{an object using \code{psVCSignal}.}

\item{...}{other parameters.}

\item{X_pred}{a matrix of \code{q} arbitrary signal vectors
of dimension \code{q} by \code{p1} for desired prediction.}

\item{t_pred}{a \code{q} vector for the varying index variable associated with \code{X_pred}.}

\item{type}{the mean value \code{type = "mu"} (default) or linear predictor
\code{type = "eta"}.}
}
\value{
\item{pred}{the estimated mean (inverse link function) (default)
or the linear predictor prediction with \code{type =
"eta"}, at signals in matrix \code{X_pred} and covariates in vector \code{t_pred}.}
}
\description{
Prediction function which returns both linear
predictor and inverse link predictions for an arbitrary matrix of
signals with their vector of companion indexing covariates (using
\code{psVCSignal} with \code{class psvcsignal}).
}
\examples{
library(fds)
data(nirc)
iindex <- nirc$x
X <- nirc$y
sel <- 50:650 # 1200 <= x & x<= 2400
X <- X[sel, ]
iindex <- iindex[sel]
dX <- diff(X)
diindex <- iindex[-1]
y <- as.vector(labc[1, 1:40]) # percent fat
t_var <- as.vector(labc[4, 1:40]) # percent flour
oout <- 23
dX <- t(dX[, -oout])
y <- y[-oout]
t_var = t_var[-oout]
Pars = rbind(c(min(diindex), max(diindex), 25, 3, 1e-7, 2),
c(min(t_var), max(t_var), 20, 3, 0.0001, 2))
fit1 <- psVCSignal(y, dX, diindex, t_var, Pars = Pars,
family = "gaussian", link = "identity", int = TRUE)
predict(fit1, X_pred = dX[1:5,], t_pred = t_var[1:5])
}
\references{
Eilers, P. H. C. and Marx, B. D. (2003). Multivariate calibration with temperature
interaction using two-dimensional penalized signal regression. \emph{Chemometrics
and Intellegent Laboratory Systems}, 66, 159–174.

Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
P-splines.} Cambridge University Press.
}
\author{
Paul Eilers and Brian Marx
}
