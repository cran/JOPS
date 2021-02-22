#' Varying-coefficient penalized signal regression using P-splines.
#'
#' @description \code{psVCSignal} is used to regress a (glm) response onto a
#'  signal such that the signal coefficients can vary over another covariate \code{t}.
#'  Anisotripic penalization of tensor product B-splines produces a 2D coefficient surface that
#'  can be sliced at \code{t}.
#'
#'  @details Support functions needed: \code{pspline_fitter}, \code{pspline_2dchecker},
#'  and \code{bbase}.
#'
#'  @import stats
#'
#' @param y a glm response vector of length \code{m}, usually continuous, binary/bimomial or counts.
#' @param X a \code{m} by \code{p1} Signal matrix of regressors.
#' @param x_index \code{p1}-vector for index of Signal (e.g. wavelength).
#' @param t_var \code{p2}-vector with other (indexing) variable in coefficient surface (e.g. temperature, depth, time).
#' @param Pars a matrix with 2 rows, each with P-spline parameters:
#' \code{min max nseg bdeg lambda pord}, for row and columns of tensor product surface; defaults are min and
#' max for \code{x_index} and \code{t_var} (resp.), \code{nseg} = 10, \code{bdeg} =3,
#' \code{lambda} = 1, \code{pord} = 2.
#' @param family the response distribution, e.g.
#' \code{"gaussian", "binomial", "poisson", "Gamma"} distribution; quotes are needed (default \code{"gaussian"}.
#' @param link the link function, one of \code{"identity"}, \code{"log"}, \code{"sqrt"},
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, \code{"reciprocal"});
#' quotes are needed (default \code{"identity"}.
#' @param m_binomial a vector of binomial trials having \code{length(y)}. Default is 1 vector for \code{family = "binomial"}, NULL otherwise.
#' @param r_gamma a vector of gamma shape parameters. Default is 1 vector for \code{family = "Gamma"}, NULL otherwise.
#' @param X_pred a matrix of signals with \code{ncol(X)} columns for prediction, default is \code{X}.
#' @param t_pred a vector for the VC indexing variable with length \code{nrow(X_pred)}, default is \code{t_var}.
#' @param y_predicted a vector for the responses associated with \code{X_pred}
#' with length \code{nrow(X_pred)} useful for CV when \code{family = "binomial"}, default is \code{NULL}.
#' @param wts a \code{m} vector of weights (default 1).
#' @param ridge_adj a small ridge penalty tuning parameter to regularize estimation (default \code{1e-8}).
#' @param int intercept set to TRUE or FALSE for intercept term.

#' @return
#' \item{pcoef}{a vector of length \code{(Pars[1,3]+Pars[1,4])}*\code{(Pars[2,3]+Pars[2,4])}
#' of estimated P-spline coefficients for tensor surface.}
#' \item{summary_predicted}{inverse link prediction vectors, and twice se bands.}
#' \item{dev}{ the deviance of fit.}
#' \item{eff_dim}{the approximate effective dimension of fit.}
#' \item{family}{the family of the response.}
#' \item{link}{the link function.}
#' \item{aic}{AIC.}
#' \item{df_resid}{approximate df residual.}
#' \item{cv}{leave-one-out standard error prediction when \code{family = "gaussian"}, NULL otherwise.}
#' \item{cv_predicted}{standard error prediction for \code{y_predict} when \code{family = "gaussian"}, NULL otherwise.}
#' \item{Pars}{design and tuning parameters; see arguments above.}
#' \item{dispersion_parm}{estimate of dispersion, \code{Dev/df_resid}.}
#' \item{summary_predicted}{inverse link prediction vectors, and twice se bands.}
#' \item{eta_predicted}{estimated linear predictor of \code{length(y)}.}
#' \item{press_mu}{leave-one-out prediction of mean when \code{family = "gaussian"}, NULL otherwise.}
#' \item{bin_percent_correct}{percent correct classification based on 0.5 cut-off when \code{family = "binomial"}, NULL otherwise.}
#' \item{Bx}{B-spline basis matrix of dimension \code{p1} by \code{n1}, along \code{x_index}.}
#' \item{By}{B-spline basis matrix of dimension \code{p2} by \code{n2}, along \code{t_var}.}
#' \item{Q}{Modified tensor basis (\code{m} by \code{(n1*n2)}) for VC signal regression.}
#' \item{yint}{the estimated y-intercept (when \code{int = TRUE}.)}
#' \item{int}{a logical variable related to use of y-intercept in model.}

#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (2003). Multivariate calibration with temperature
#' interaction using two-dimensional penalized signal regression. \emph{Chemometrics
#' and Intellegent Laboratory Systems}, 66, 159–174.
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

"psVCSignal" <-
  function(y, X, x_index, t_var, Pars= rbind(c(min(x_index),
             max(x_index), 10, 3, 1, 2), c(min(t_var),
             max(t_var), 10, 3, 1, 2)),
             family = "gaussian", link = "default",
             m_binomial = 1 + 0*y, wts = 1 + 0*y , r_gamma = 1 + 0*y,
             X_pred = X, t_pred = t_var, y_predicted = NULL,
             ridge_adj = 1e-8, int = TRUE) {
    #' @export

    # Dimension, weight, parameter defaults
    m <- length(y)

    # 	Check to see if any argument is improperly defined
    parms <- pspline2d_checker(family, link, Pars[1, 4], Pars[2, 4], Pars[1, 6], Pars[2, 6],
                               Pars[1, 3], Pars[2, 3], Pars[1, 5], Pars[2, 5], ridge_adj, wts)
    family <- parms$family
    link <- parms$link
    ridge_adj <- parms$ridge_adj
    wts <- parms$wts
    Pars[1, 3:6] <- c(parms$nseg1, parms$bdeg1, parms$lambda1, parms$pord1)
    Pars[2, 3:6] <- c(parms$nseg2, parms$bdeg2, parms$lambda2, parms$pord2)

    # Prepare bases for estimation
    X <- as.matrix(X)
    Bx <- bbase(x_index, Pars[1, 1 ], Pars[1, 2], Pars[1, 3], Pars[1, 4])
    U <- X %*% Bx
    By <- bbase(t_var, Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
    n1 <- ncol(Bx)
    n2 <- ncol(By)

    # Compute tensor products
    XB1 <- kronecker(U, t(rep(1, n2)))

    B1 <- kronecker(Bx, t(rep(1, n2)))
    B2 <- kronecker(t(rep(1, n1)), By)

    # Modified VC tensor product basis
    Q <- XB1 * B2
    Qnoint <- Q

    # Prepare penalty matrices
    d1 <- Pars[1, 6]
    D1 <- diff(diag(n1), diff = d1)
    lambda1 <- Pars[1, 5]
    P1 <- sqrt(lambda1) * kronecker(D1, diag(n2))
    d2 <- Pars[2, 6]
    D2 <- diff(diag(n2), diff = d2)
    lambda2 <- Pars[2, 5]
    P2 <- sqrt(lambda2) * kronecker(diag(n1), D2)
    Pen <- rbind(P1, P2)


    # Prepare ridge penalty
    p_ridge <- NULL
    if (ridge_adj > 0) {
      nix_ridge <- rep(0, n1 * n2)
      p_ridge <- sqrt(ridge_adj) * diag(n1 * n2)
    }

    # Data augmentation and regression, with or without intercept
    z1 <- rep(0, n2 * (n1 - d1))
    z2 <- rep(0, n1 * (n2 - d2))
    n_col <- ncol(Q)
    if (int) {
      Q <- cbind(rep(1, nrow(Q)), Q)
      Pen <- cbind(rep(0, nrow(Pen)), Pen)
      if (ridge_adj > 0) {
        p_ridge <- cbind(rep(0, nrow(p_ridge)), p_ridge)
      }
    }

    # Apply the penalized GLM fitter on for y on Q
    ps_fit <- pspline_fitter(y, B = Q, family, link, P = Pen, P_ridge = p_ridge,
                             wts, m_binomial, r_gamma)
    mu <- ps_fit$mu
    pcoef <- ps_fit$coef
    bin_percent_correct <- NULL

    # Percent correct classification for binomial
    if (family == "binomial") {
      pcount <- 0
      p_hat <- mu / m_binomial
      for (ii in 1:m) {
        if (p_hat[ii] > 0.5) {
          count <- y[ii]
        }
        if (p_hat[ii] <= 0.5) {
          count <- m_binomial[ii] - y[ii]
        }
        count <- pcount + count
        pcount <- count
      }
      bin_percent_correct <- count / sum(m_binomial)
    }

    # GLM weights
    w <- ps_fit$w

    # Hat diagonal for ED
    h <- hat(ps_fit$f$qr, intercept = FALSE)[1:m]
    trace <- eff_dim <- sum(h)

    # Deviance, by family
    dev_ = dev_calc(family, y, mu, m_binomial, r_gamma)
    dev = dev_$dev
    dispersion_parm = dev_$dispersion_parm

    # Leave-one-out CV for Gaussian	response
    cv <- press_mu <- press_e <- NULL
    if (family == "gaussian") {
      dev <- sum(ps_fit$f$residuals[1:m]^2)
      dispersion_parm <- dev / (m - trace)
      press_e <- ps_fit$f$residuals[1:m] / (1 - h)
      cv <- sqrt(sum((press_e)^2) / (m))
      press_mu <- y - press_e
    }

    # AIC
    aic <- dev + 2 * trace

    # Weight augmentation
    w_aug <- c(w, (c(z1, z2) + 1))

    # Isolating intercept
    if (int) {
      yint <- ps_fit$coef[1]
    }
    if (!int) {
      yint <- NULL
    }

    # External prediction
    summary_predicted <- NULL
    cv_predicted <- eta_predicted <- avediff_pred <- NULL
    if (!missing(t_pred)) {
      q <- length(t_pred)
      Up <- as.matrix(X_pred) %*% Bx
      Byp <- bbase(t_pred, Pars[2, 1 ], Pars[2, 2], Pars[2, 3], Pars[2, 4])
      B1p <- kronecker(Up, t(rep(1, n2)))
      B2p <- kronecker(t(rep(1, n1)), Byp)
      Qp <- B1p * B2p
      if (!int) {
        eta_predicted <- Qp %*% pcoef
      }
      if (int) {
        one_xpred_b <- cbind(rep(1, q), Qp)
        eta_predicted <- Qp %*% pcoef[2:length(pcoef)] + yint
      }

      summary_predicted <- eta_predicted
      if (!missing(y_predicted)) {
        if (family == "gaussian") {
          cv_predicted <- sqrt(sum((y_predicted - eta_predicted)^2) / (length(
            y_predicted
          )))
          avediff_pred <- (sum(y_predicted - eta_predicted)) / length(y_predicted)
        }
      }

      # Future prediction, percent correct classification
      bin_percent_correct <- NULL
      if (link == "logit") {
        summary_predicted <- 1 / (1 + exp(-summary_predicted))
        pcount <- 0
        p_hat <- exp(eta_predicted) / (1 + exp(eta_predicted))
        if (!missing(y_predicted)) {
          for (ii in 1:length(eta_predicted)) {
            if (p_hat[ii] > 0.5) {
              count <- y_predicted[ii]
            }
            if (p_hat[ii] <= 0.5) {
              count <- 1 - y_predicted[ii]
            }
            count <- pcount + count
            pcount <- count
          }
          bin_percent_correct <- count / length(y_predicted)
        }
      }

      # Future prediction with inverse link
      summary_predicted <- inverse_link(summary_predicted, link)
    }

    # Return list
    P <- list(
      pcoef = pcoef, Pars = Pars, cv = cv, eff_dim = eff_dim, yint = yint, int = int,
      bin_percent_correct = bin_percent_correct, family = family, link = link, aic = aic, dev = dev, df_resid
      = m - trace, dispersion_parm = dispersion_parm, mu = mu, press_mu = press_mu,
      summary_predicted = summary_predicted, cv_predicted = cv_predicted, eta_predicted =
        eta_predicted, avediff_pred = avediff_pred, ridge_adj = ridge_adj, Q = Q, Bx = Bx, By = By)
    class(P) <- "psvcsignal"
    return(P)
  }
