#' Predict function for \code{psNormal}, \code{psBinomial}, \code{psPoisson}
#'
#' @description Prediction function which returns both linear
#' predictor and inverse link predictions at arbitrary data locations
#' (using \code{psNormal}, \code{psBinomial}, \code{psPoisson} with \code{class pspfit}).
#'
#' @param object an object using \code{psNormal}, \code{psBinomial}, or \code{psPoisson} .
#' @param ... other parameters.
#' @param x a scalar or vector of arbitrary \code{x} locations for
#' desired prediction.
#' @param type the mean value \code{type = "mu"} (default) or linear predictor
#' \code{type = "eta"}.

#' @return
#' \item{pred}{the estimated mean (inverse link function) (default)
#' or the linear predictor prediction with \code{type =
#' "eta"}, at arbitary \code{x} locations.}

#' @author Paul Eilers and Brian Marx
#' @references  Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' library(JOPS)
#' library(boot)
#'
#' # Extract the data
#' Count <- hist(coal$date, breaks = c(1851:1963), plot = FALSE)$counts
#' Year <- c(1851:1962)
#' xl <- min(Year)
#' xr <- max(Year)
#'
#' # Poisson smoothing
#' nseg <- 20
#' bdeg <- 3
#' fit1 <- psPoisson(Year, Count, xl, xr, nseg, bdeg, pord = 2, lambda = 1)
#' names(fit1)
#' plot(fit1, xlab = "Year", ylab = "Count", se = 2)
#' predict(fit1, x = fit1$x[1:5])
#' predict(fit1, x = fit1$x[1:5], type = "eta")
#' @export
predict.pspfit = function(object,..., x, type = "mu"){

  # Retrieve design parameters
  ps = object
  nseg = ps$nseg
  bdeg = ps$bdeg
  xl = ps$xl
  xr = ps$xr
  pcoef = ps$pcoef

  # Code block for psNormal
  if (ps$family == 'gaussian') {
    Bnew = bbase(x, xl, xr, nseg, bdeg)
    pred <- as.vector(Bnew%*%pcoef)
    return(pred)
  }

  if (ps$family == 'binomial') {
    Bnew = bbase(x, xl, xr, nseg, bdeg)
    pred <- as.vector(Bnew%*%pcoef)
    if(type == 'mu') {
      pred = 1/(1 + exp(-pred))
    }
    return(pred)

  }

  if (ps$family == 'poisson') {
    Bnew = bbase(x, xl, xr, nseg, bdeg)
    pred <- as.vector(Bnew%*%pcoef)
    if(type == 'mu'){
      pred = exp(pred)
    }
    return(pred)
  }
}

#' Predict function for \code{psSignal}
#'
#' @description Prediction function which returns both linear
#' predictor and inverse link predictions, for an arbitrary matrix of signals
#' (using \code{psSignal} with \code{class pssignal}).
#'
#' @param object an object using \code{psSignal}.
#' @param ... other parameters.
#' @param X_pred a matrix of arbitrary signals with \code{ncol(X) == length(x_index)} locations for
#' desired prediction.
#' @param type the mean value \code{type = "mu"} (default) or linear predictor
#' \code{type = "eta"}.

#' @return
#' \item{pred}{the estimated mean (inverse link function) (default)
#' or the linear predictor prediction with \code{type =
#' "eta"}, for a matrix of signals in \code{X_pred}.}

#' @author Paul Eilers and Brian Marx
#' @references  Marx, B.D. and Eilers, P.H.C. (1999). Generalized linear regression for sampled signals and
#'        curves: A P-spline approach. \emph{Technometrics}, 41(1): 1-13.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' library(JOPS)
#' # Get the data
#' library(fds)
#' data(nirc)
#' iindex=nirc$x
#' X=nirc$y
#' sel= 50:650 #1200 <= x & x<= 2400
#' X=X[sel,]
#' iindex=iindex[sel]
#' dX=diff(X)
#' diindex=iindex[-1]
#' y=as.vector(labc[1,1:40])
#' oout=23
#' dX=t(dX[,-oout])
#' y=y[-oout]
#' fit1 = psSignal(y, dX, diindex, nseg = 25,lambda = 0.0001)
#' predict(fit1, X_pred = dX[1:5, ])
#' predict(fit1, X_pred = dX[1:5, ], type = 'eta')
#' @export
predict.pssignal = function(object,..., X_pred, type = "mu"){
  # Retrieve design parameters
  pssig = object
  B = pssig$B
  pcoef = pssig$coef
  yint = pssig$y_int
  link = pssig$link
  family = pssig$family
  X = X_pred

  pred <- as.vector((X%*%B%*%pcoef[-1]) + yint)
  if(type == 'mu'){
    pred = inverse_link(pred, link)
  }
  return(pred)
}

#' Predict function for \code{ps2DNormal}
#'
#' @description Prediction function which returns linear
#'  predictions at arbitrary (x, y) data locations (using \code{ps2DNormal}
#'  with \code{class ps2dnormal}).
#'
#' @param object an object using ps2DNormal.
#' @param ... other parameters.
#' @param XY a matrix of arbitrary (\code{x}, \code{y}) locations for
#' desired prediction.
#'
#' @return
#' \item{pred}{the estimated mean at (\code{x}, \code{y}) locations, in \code{XY}.}

#' @author Paul Eilers and Brian Marx
#' @references  Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @examples
#' library(SemiPar)
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
#' xpars <- c(xlo, xhi, 10, 3, 0.01, 1)
#' ypars <- c(ylo, yhi, 10, 3, 0.1, 1)
#' Pars1 <- rbind(xpars, ypars)
#' fit <- ps2DNormal(cbind(x, y, z), Pars = Pars1)
#' predict(fit, XY = cbind(x, y)[1:5, ])
#' @export
predict.ps2dnormal = function(object,..., XY){

  # Retrieve design parameters
  ps2dnor = object
  Pars = ps2dnor$Pars
  p1 <- Pars[1, ]
  p2 <- Pars[2, ]
  nx = p1[3] + p1[4]
  ny = p2[3] + p2[4]
  pcoef = ps2dnor$coef

  # Prediction
  Bxp <- bbase(as.vector(XY[, 1]), p1[1], p1[2], p1[3], p1[4])
  Byp <- bbase(as.vector(XY[, 2]), p2[1], p2[2], p2[3], p2[4])
  B1p <- kronecker(Bxp, t(rep(1, ny)))
  B2p <- kronecker(t(rep(1, nx)), Byp)
  Bp <- B1p * B2p
  pred <- as.vector(Bp %*% pcoef)
  return(pred)
}

#' Predict function for \code{ps2DSignal}
#'
#' @description Prediction function which returns both linear
#' predictor and inverse link predictions for arbitrary 2D signals (using
#' \code{ps2DSignal} with \code{class ps2dsignal}).
#'
#' @param object an object using \code{ps2DSignal}.
#' @param ... other parameters.
#' @param M_pred a matrix of \code{q} arbitrary "stacked" or "unfolded" signal matrices
#' of dimension (\code{q} by \code{p1}) by \code{p2} or \code{q} by (\code{p1}
#' by \code{p2}, respectively,
#' for desired prediction (default "unfolded").
#' @param M_type "stacked" or "unfolded" (default).
#' @param type the mean value \code{type = "mu"} (default) or linear predictor
#' \code{type = "eta"}.

#' @return
#' \item{pred}{the estimated mean (inverse link function)
#' or the linear predictor prediction with \code{type =
#' "eta"}, for arbitary 2D signals in \code{M_pred}.}

#' @author Paul Eilers and Brian Marx

#' @references Marx, B.D. and Eilers, P.H.C. (2005).
#' Multidimensional penalized signal regression, \emph{Technometrics}, 47: 13-22.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' library(fields)
#' library(JOPS)
#'
#' # Get the data
#' x0 <- Sugar$X
#' x0 <- x0 - apply(x0, 1, mean) # center Signal
#' y <- as.vector(Sugar$y[, 3]) # Response is Ash
#'
#' # Inputs for two-dimensional signal regression
#' nseg <- c(7, 37)
#' pord <- c(3, 3)
#' min_ <- c(230, 275)
#' max_ <- c(340, 560)
#' M1_index <- rev(c(340, 325, 305, 290, 255, 240, 230))
#' M2_index <- seq(from = 275, to = 560, by = .5)
#' p1 <- length(M1_index)
#' p2 <- length(M2_index)
#'
#' # Fit optimal model based on LOOCV
#' opt_lam <- c(8858.6679, 428.1332) # Found via svcm
#' Pars_opt <- rbind(
#'   c(min_[1], max_[1], nseg[1], 3, opt_lam[1], pord[1]),
#'   c(min_[2], max_[2], nseg[2], 3, opt_lam[2], pord[2])
#' )
#' fit <- ps2DSignal(y, x0, p1, p2, "unfolded", M1_index, M2_index,
#'   Pars_opt,int = TRUE, ridge_adj = 0.0001,
#'   M_pred = x0 )
#'
#'
#' predict(fit, M_pred= x0, type = "mu", M_type = "unfolded")

#' @export
predict.ps2dsignal = function(object,..., M_pred, M_type = "unfolded",
                              type = "mu")
  {

# Retrieve design parameters
  ps2dsig = object
  M1index = ps2dsig$M1index
  M2index = ps2dsig$M2index
  p1 = length(M1index)
  p2 = length(M2index)
  Ahat = ps2dsig$Ahat
  int = ps2dsig$int
  yint = ps2dsig$yint
  link = ps2dsig$link
# External prediction
  if (nrow(M_pred)>1) {
    q <- nrow(M_pred)
    X_p <- as.matrix(M_pred)
    if (M_type == "stacked") {
      M_pred <- as.matrix(M_pred)
      q <- nrow(M_pred) / p1
      x <- as.vector(t(M_pred))
      X_p <- matrix(x, q, p1 * p2, byrow = TRUE)
    }
      pred <- X_p %*% as.vector(Ahat)

    if (int) {
      pred <- X_p %*% Ahat + yint
      }
  }
  if(type == 'mu'){
    pred = inverse_link(pred, link)
  }
return(pred)
}

#' Predict function for \code{psVCSignal}
#'
#' @description Prediction function which returns both linear
#' predictor and inverse link predictions for an arbitrary matrix of
#' signals with their vector of companion indexing covariates (using
#' \code{psVCSignal} with \code{class psvcsignal}).
#'
#' @param object an object using \code{psVCSignal}.
#' @param ... other parameters.
#' @param X_pred a matrix of \code{q} arbitrary signal vectors
#' of dimension \code{q} by \code{p1} for desired prediction.
#' @param t_pred a \code{q} vector for the varying index variable associated with \code{X_pred}.
#' @param type the mean value \code{type = "mu"} (default) or linear predictor
#' \code{type = "eta"}.

#' @return
#' \item{pred}{the estimated mean (inverse link function) (default)
#' or the linear predictor prediction with \code{type =
#' "eta"}, at signals in matrix \code{X_pred} and covariates in vector \code{t_pred}.}

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
#' predict(fit1, X_pred = dX[1:5,], t_pred = t_var[1:5])

#' @export
predict.psvcsignal = function(object,..., X_pred, t_pred, type = "mu")
{
  # Retrieve design parameters
  psvcsig = object
  M1index = psvcsig$M1index
  M2index = psvcsig$M2index
  p1 = length(M1index)
  p2 = length(M2index)
  pcoef = psvcsig$pcoef
  int = psvcsig$int
  yint = psvcsig$yint
  link = psvcsig$link
  Pars = psvcsig$Pars
  Bx = psvcsig$Bx
  n1 = ncol(Bx)
  By = psvcsig$By
  n2 = ncol(By)


  # External prediction
  q <- length(t_pred)
  Up <- as.matrix(X_pred) %*% Bx
  Byp <- bbase(t_pred, Pars[2, 1 ], Pars[2, 2], Pars[2, 3], Pars[2, 4])
  B1p <- kronecker(Up, t(rep(1, n2)))
  B2p <- kronecker(t(rep(1, n1)), Byp)
  Qp <- B1p * B2p

  if (!int) {
    pred <- Qp %*% pcoef
  }
  if (int) {
    one_xpred_b <- cbind(rep(1, q), Qp)
    pred <- Qp %*% pcoef[2:length(pcoef)] + yint
  }
  if(type == 'mu'){
    pred = inverse_link(pred, link)
  }
  return(pred)
}


#' Predict function for \code{ps2DGLM}
#'
#' @description Prediction function which returns both linear
#' predictor and inverse link predictions at arbitrary (x, y) data locations
#' (using \code{ps2DGLM} with \code{class ps2dglm}).
#'
#' @param object an object using \code{ps2DGLM}.
#' @param ... other parameters.
#' @param XY a matrix of arbitrary (\code{x}, \code{y}) locations for
#' desired prediction.
#' @param type the mean value \code{type = "mu"} (default) or linear predictor
#' \code{type = "eta"}.
#'
#' @return
#' \item{pred}{the estimated mean (inverse link function) (default)
#' or the linear predictor prediction with \code{type =
#' "eta"}, for arbitary (x, y) locations in \code{XY}.}

#' @author Paul Eilers and Brian Marx

#' @references  Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
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
#'   family = "binomial", link = "logit")
#' predict(fit, XY = cbind(Start, Age)[1:5,])


#' @export
predict.ps2dglm = function(object,..., XY, type = "mu")
{

  # Retrieve design parameters
  ps2dg = object
  Pars = ps2dg$Pars
  pcoef = ps2dg$pcoef
  link = ps2dg$link
  n1 = Pars[1, 3] + Pars[1, 4]
  n2 = Pars[2, 3] + Pars[2, 4]
  # External prediction
  if(is.vector(XY)) {XY = t(XY)}
  Bxp <- bbase(XY[, 1], Pars[1, 1], Pars[1, 2], Pars[1, 3], Pars[1, 4])
  Byp <- bbase(XY[, 2], Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
  B1p <- kronecker(Bxp, t(rep(1, n2)))
  B2p <- kronecker(t(rep(1, n1)), Byp)
  Bp <- B1p * B2p
  pred <- Bp %*% pcoef

  if(type == 'mu'){
    pred = inverse_link(pred, link)
  }
  return(pred)
}

#' Predict function for \code{sim_psr}
#'
#' @description Prediction function which returns single-index
#' inverse link linear predictions at arbitrary data locations (using
#' \code{sim_psr} with \code{class simpsr}).
#'
#' @param object an object using \code{sim_psr}.
#' @param ... other parameters.
#' @param X_pred a matrix of arbitrary signals with \code{ncol(X_pred) = length(x_index)} locations for
#' desired prediction.
#'
#' @return
#' \item{pred}{the estimated (inverse single-index) mean for the signals in \code{X_pred}.}

#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C., B. Li, B.D. Marx (2009).
#' Multivariate calibration with single-index signal regression,
#' \emph{Chemometrics and Intellegent Laboratory Systems}, 96(2), 196-202.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @examples
#' library(JOPS)
#' # Get the data
#' library(fds)
#' data(nirc)
#' iindex <- nirc$x
#' X <- nirc$y
#' sel <- 50:650 # 1200 <= x & x<= 2400
#' X <- X[sel, ]
#' iindex <- iindex[sel]
#' dX <- diff(X)
#' diindex <- iindex[-1]
#' y <- as.vector(labc[1, 1:40])
#' oout <- 23
#' dX <- t(dX[, -oout])
#' y <- y[-oout]
#'
#' pords <- c(2, 2)
#' nsegs <- c(27, 7)
#' bdegs = c(3, 3)
#' lambdas <- c(1e-6, .1)
#' max_iter <- 100
#'
#' # Single-index model
#' fit <- sim_psr(y, dX, diindex, nsegs, bdegs, lambdas, pords,
#'              max_iter)
#' predict(fit, X_pred = dX)
#'
#' @export
#'
predict.simpsr <- function(object,..., X_pred) {
  sim = object
  B = sim$B
  yint = sim$yint
  alpha = sim$alpha
  f_fit = sim$f_fit
  eta <- as.vector(X_pred %*% B %*% alpha) + yint
  pred <- predict(f_fit, x = eta)
  return(pred)
}

#' Predict function for \code{sim_vcpsr}
#'
#' @description Prediction function which returns varying-coefficient single-index
#' inverse link linear predictions at arbitrary data locations (using \code{sim_vcpsr} with
#' \code{class simvcpsr}).
#'
#' @param object an object using \code{sim_vcpsr}.
#' @param ... other parameters.
#' @param X_pred a matrix of arbitrary signals with \code{ncol(X_pred) = length(x_index)} locations for
#' desired prediction.
#' @param t_pred a \code{q} vector for the VC index variable associated with \code{X_pred}.
#' @return
#' \item{pred}{the estimated (inverse single-index) mean for the signals in the matrix \code{X_pred},
#' with the companion vector of indexing covariates in \code{t_pred}.}

#' @author Paul Eilers and Brian Marx
#' @references Marx, B. D. (2015). Varying-coefficient single-index signal
#' regression. \emph{Chemometrics and Intellegent Laboratory Systems}, 143, 111–121.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
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
#' predict(fit, X_pred = x, t_pred = t_var)


#' @export
#'
predict.simvcpsr <- function(object,..., X_pred, t_pred) {
  simvc = object
  Pars = simvc$Pars
  Bx = simvc$Bx
  yint = simvc$yint
  alpha = simvc$alpha
  By = simvc$By
  fit2D = simvc$fit2D

  U <- X_pred %*% Bx
  By <- bbase(t_pred, Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
  n1 <- ncol(Bx)
  n2 <- ncol(By)

  # Compute tensor products
  XB1 <- kronecker(U, t(rep(1, n2)))

  B1 <- kronecker(Bx, t(rep(1, n2)))
  B2 <- kronecker(t(rep(1, n1)), By)

  # Modified VC tensor product basis
  Q <- XB1 * B2

  eta <- as.vector(Q %*% alpha) + yint
  pred <- predict(fit2D, XY = cbind(eta, t_pred))
  return(pred)
}
