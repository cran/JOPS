#' Smooth signal (multivariate calibration) regression using P-splines.
#' @description Smooth signal (multivariate calibration) regression using P-splines.
#'
#' @details Support functions needed: \code{pspline_fitter}, \code{bbase} and \code{pspline_checker}.
#'
#' @import stats
#'
#' @author Brian Marx
#'
#' @param y a (glm) response vector, usually continuous, binomial or count data.
#' @param x_signal  a matrix of continuous regressor with \code{nrow(x_signal) == length(y)}, often
#' a discrete digitization of a signal or histogram or time series.
#' @param x_index a vector to of length \code{ncol(x_signal) == p}, associated with the
#' ordering index of the signal. Default is \code{1:ncol(x_signal)}.
#' @param family the response distribution, e.g.
#' \code{"gaussian", "binomial", "poisson", "Gamma"} distribution; quotes are needed. Default is \code{"gaussian"}.
#' @param link the link function, one of \code{"identity"}, \code{"log"}, \code{"sqrt"},
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, \code{"reciprocal"};
#' quotes are needed (default \code{"identity"}).
#' @param nseg the number of evenly spaced segments between \code{xl} and \code{xr} (default 10).
#' @param bdeg the degree of the basis, usually 1, 2, or 3 (defalult).
#' @param lambda the (positive) tuning parameter for the penalty (default 1).
#' @param pord the order of the difference penalty, usually 1, 2, or 3 (defalult).
#' @param m_binomial a vector of binomial trials having length(y); default is 1 vector for \code{family = "binomial"}, NULL otherwise.
#' @param r_gamma a vector of gamma shape parameters. Default is 1 vector for \code{family = "Gamma"}, NULL otherwise.
#' @param wts the weight vector of \code{length(y)}; default is 1.
#' @param x_predicted a matrix of external signals to yield external prediction.
#' @param y_predicted a vector of responses associated
#' with \code{x_predicted} which are used to calculate standard error of external prediction. Default is NULL.
#' @param ridge_adj A ridge penalty tuning parameter, which can be set to small value, e.g. \code{1e-8} to stabilize estimation, (default 0).
#' @param int set to TRUE or FALSE to include intercept term in linear predictor (default TRUE).

#' @return
#' \item{coef}{a vector with \code{length(n)} of estimated P-spline coefficients.}
#' \item{mu}{a vector with \code{length(m)} of estimated means.}
#' \item{eta}{a vector of \code{length(m)} of estimated linear predictors.}
#' \item{B}{the B-spline basis (for the coefficients), with dimension \code{p} by \code{n}.}
#' \item{deviance}{the deviance of fit.}
#' \item{eff_df}{the approximate effective dimension of fit.}
#' \item{aic}{AIC.}
#' \item{df_resid}{approximate df residual.}
#' \item{beta}{a vector of length \code{p}, containing estimated smooth signal coefficients.}
#' \item{std_beta}{a vector of length \code{p}, containing standard errors of smooth signal coefficients.}
#' \item{cv}{leave-one-out standard error prediction, when \code{family = "gaussian"}.}
#' \item{cv_predicted}{standard error prediction for \code{y_predict}, when \code{family = "gaussian"}, NULL otherwise.}
#' \item{nseg}{the number of evenly spaced B-spline segments.}
#' \item{bdeg}{the degree of B-splines.}
#' \item{pord}{the order of the difference penalty.}
#' \item{lambda}{the positive tuning parameter.}
#' \item{family}{the family of the response.}
#' \item{link}{the link function.}
#' \item{y_intercept}{the estimated y-intercept (when \code{int = TRUE}.)}
#' \item{int}{a logical variable related to use of y-intercept in model.}
#' \item{dispersion_param}{estimate of dispersion, \code{Dev/df_resid}.}
#' \item{summary_predicted}{inverse link prediction vectors, and twice se bands.}
#' \item{eta_predicted}{estimated linear predictor of \code{length(y)}.}
#' \item{press_mu}{leave-one-out prediction of mean, when \code{family = "gaussian"}, NULL otherwise.}
#' \item{bin_percent_correct}{percent correct classification based on 0.5 cut-off, when \code{family = binomial}, NULL otherwise.}
#' \item{x_index}{a vector to of length \code{ncol(x_signal) == p}, associated with the ordering of the signal.}

#' @references  Marx, B.D. and Eilers, P.H.C. (1999). Generalized linear regression for sampled signals and
#' curves: A P-spline approach. \emph{Technometrics}, 41(1): 1-13.
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
#' y <- as.vector(labc[1, 1:40]) # percent fat
#' oout <- 23
#' dX <- t(dX[, -oout])
#' y <- y[-oout]
#' fit1 <- psSignal(y, dX, diindex, nseg = 25, bdeg = 3, lambda = 0.0001,
#' pord = 2, family = "gaussian", link = "identity", x_predicted = dX, int = TRUE)
#' plot(fit1, xlab = "Coefficient Index", ylab = "ps Smooth Coeff")
#' title(main = "25 B-spline segments with tuning = 0.0001")
#' names(fit1)
#' @export
"psSignal" <-
  function(y, x_signal, x_index = c(1: ncol(x_signal)), nseg = 10,
           bdeg = 3, pord = 3, lambda = 1, wts = 1+ 0*y, family =
           "gaussian", link = "default", m_binomial = 1 + 0*y,
           r_gamma = wts, y_predicted = NULL, x_predicted
             = x_signal, ridge_adj = 0, int = TRUE) {
    x <- x_index
    n <- length(y)

    # 	Check to see if any argument is improperly defined
    parms <- pspline_checker(
      family, link, bdeg, pord,
      nseg, lambda, ridge_adj, wts
    )
    family <- parms$family
    link <- parms$link
    q <- parms$bdeg
    d <- parms$pord
    ridge_adj <- parms$ridge_adj
    lambda <- parms$lambda
    nseg <- parms$nseg
    wts <- parms$wts

    # Construct basis
    xl <- min(x)
    xr <- max(x)
    xmax <- xr + 0.01 * (xr - xl)
    xmin <- xl - 0.01 * (xr - xl)
    dx <- (xmax - xmin) / nseg
    knots <- seq(xmin - q * dx, xmax + q * dx, by = dx)
    b <- bbase(x, xl, xr, nseg, q)
    n_col <- ncol(b)

    # Define the ridge and difference penalty matrices
    p_ridge <- NULL
    if (ridge_adj > 0) {
      nix_ridge <- rep(0, n_col)
      p_ridge <- sqrt(ridge_adj) * diag(rep(1, n_col))
    }
    p <- diag(n_col)
    p <- diff(p, diff = d)

    # Set up data augmentation
    p <- sqrt(lambda) * p
    nix <- rep(0, n_col - d)

    # Construct effective regresssors U=XB, with intercept
    x_signal <- as.matrix(x_signal)
    xb <- x_signal %*% as.matrix(b)
    if (int) {
      xb <- cbind(rep(1, n), xb)
      p <- cbind(rep(0, nrow(p)), p)
      if (ridge_adj > 0) {
        p_ridge <- cbind(rep(0, nrow(p_ridge)), p_ridge)
      }
    }

    # Apply the penalized GLM fitter on for y on U=XB
    ps_fit <- pspline_fitter(y, B = xb, family, link, P = p,
      P_ridge = p_ridge, wts, m_binomial, r_gamma )
    mu <- ps_fit$mu
    coef <- ps_fit$coef
    eta <- ps_fit$eta

    # Percent classification for Binomial
    bin_percent_correct <- NULL
    if (family == "binomial") {
      pcount <- 0
      p_hat <- mu / m_binomial
      for (ii in 1:n) {
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
    w_aug <- c(w, (nix + 1))

    # Hat diagonal for ED
    qr = ps_fit$f$qr
    h <- hat(qr, intercept = FALSE)[1:n]
    trace <- sum(h)

    # SE bands on a grid using QR
    R <- qr.R(qr)
    bread = chol2inv(R)

    # Deviance, by family
    dev_ = dev_calc(family, y, mu, m_binomial, r_gamma)
    dev = dev_$dev
    dispersion_parm = dev_$dispersion_parm

    # Leave-one-out CV for Gaussian	response
    cv <- press_mu <- press_e <- NULL
    if (family == "gaussian") {
      dev <- sum(ps_fit$f$residuals[1:n]^2)
      dispersion_parm <- dev / (n - trace)
      press_e <- ps_fit$f$residuals[1:n] / (1 - h)
      cv <- sqrt(sum((press_e)^2) / (n))
      press_mu <- y - press_e
    }

    # AIC
    aic <- dev + 2 * trace

    # Smooth coefficient vectors
    if (int) {
      yint <- ps_fit$coef[1]
      beta <- b %*% (as.vector(ps_fit$f$coef)[2:(n_col + 1)])
    }
    if (!int) {
      yint <- NULL
      beta <- b %*% as.vector(ps_fit$f$coef)
    }

    # Standard error bands for smooth coefficient vector
    ii = c((1+int):(n_col+int))
    var_beta = b%*% bread[ii, ii]%*%t(b)
    stdev_beta = as.vector(sqrt(dispersion_parm*diag(var_beta)))

    # Future prediction with se bands
    summary_predicted <- NULL
    cv_predicted <- eta_predicted <- NULL
    if (!missing(x_predicted)) {
      x_predicted <- as.matrix(x_predicted)
      if (!int) {
        if (ncol(x_predicted) > 1) {
          eta_predicted <- x_predicted %*% beta
          var_pred <- x_predicted %*% b %*% bread %*% t(
            x_predicted %*% b
          )
        }
        if (ncol(x_predicted) == 1) {
          eta_predicted <- t(x_predicted) %*% beta
          var_pred <- t(x_predicted) %*% b %*% bread %*%
            t(b) %*% x_predicted
        }
      }
      if (int) {
         dim_xp <- nrow(x_predicted)
        if (ncol(x_predicted) > 1) {
          one_xpred_b <- cbind(rep(1, dim_xp), (x_predicted %*% b))
          eta_predicted <- x_predicted %*% beta + yint
          var_pred <- one_xpred_b %*% bread %*% t(one_xpred_b)
        }
        if (ncol(x_predicted) == 1) {
          one_xpred_b <- cbind(1, t(x_predicted) %*% b)
          eta_predicted <- t(x_predicted) %*% beta + yint
          var_pred <- (one_xpred_b) %*% bread %*% t(
            one_xpred_b
          )
        }
      }
      stdev_pred <- as.vector(sqrt(diag(var_pred)))
      stdev_pred <- sqrt(dispersion_parm) * stdev_pred
      pivot <- as.vector(2 * stdev_pred)
      upper <- eta_predicted + pivot
      lower <- eta_predicted - pivot
      summary_predicted <- cbind(lower, eta_predicted, upper)
      if (!missing(y_predicted)) {
        if (family == "gaussian") {
          cv_predicted <- sqrt(sum((y_predicted -
            eta_predicted)^2) / (length(y_predicted)))
        }
      }
      # Future prediction, percent correct classification family binomial
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

      # Future prediction and lower, upper CI for inverse link
      summary_predicted <- inverse_link(summary_predicted, link)
      if (link == "reciprocal") {
        summary_predicted <- summary_predicted[, 3:1]
      }
      summary_predicted <- as.matrix(summary_predicted)
      dimnames(summary_predicted) <- list(NULL, c(
        "-2std_Lower",
        "Predicted", "+2std_Upper"
      ))
    }

    # Output list
    llist <- list(
      b= b, B = b, coef = coef, y_intercept = yint, int = int, x_index = x_index,
      x_signal = x_signal, y = y, press_mu = press_mu, bin_percent_correct = bin_percent_correct,
      family = family, link = link, nseg = nseg, pord = d, bdeg = q,
      lambda = lambda, aic = aic, deviance = dev, eff_df = trace - 1,
      df_resid = n - trace + 1, bin_percent_correct = bin_percent_correct,
      dispersion_param = dispersion_parm, summary_predicted = summary_predicted,
      eta_predicted = eta_predicted, cv_predicted = cv_predicted, cv = cv,
      mu = mu, eta = eta, beta = beta, stdev_beta = stdev_beta
    )
    class(llist) <- "pssignal"
    return(llist)
  }
