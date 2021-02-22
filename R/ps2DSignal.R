#' Two-dimensional penalized signal regression using P-splines.
#'
#' @description \code{ps2DSignal} is a function used to regress a (glm) response onto a two-dimensional
#' signal or image, with aniosotripic penalization of tensor product P-splines.
#'
#' @details Support functions needed: \code{pspline_fitter}, \code{bbase}, and \code{pspline_2dchecker}.
#'
#' @import stats
#'
#' @param y a response vector of length \code{m}, usually continuous, binary/bimomial or counts.
#' @param p1 the row dimension of the image.
#' @param p2 the column dimension of the image.
#' @param M_type "stacked" (signal as matrix) or "unfolded" (signal as vector).
#' @param M1_index an index of length \code{p1} for rows of regressor matrix (default is a simple sequence).
#' @param M2_index an index of length \code{p2} for columns of regressor matrix (default is a simple sequence).
#' @param Pars a matrix of 2 rows, where the first and second row
#' sets the P-spline paramters for \code{x} and \code{y}, respectively.
#' Each row consists of: \code{min max nseg bdeg lambda pord}.
#' The \code{min} and \code{max} set the ranges, \code{nseg} (default 10)
#' is the number of evenly spaced segments between \code{min} and \code{max},
#' \code{bdeg} is the degree of the basis (default 3 for cubic),
#' \code{lambda} is the (positive) tuning parameter for the penalty (default 1),
#' \code{pord} is the number for the order of the difference penalty (default 2).
#' @param M The signal/image regressors, which are either "stacked" or "unfolded",
#' with dimensions (\code{m} * \code{p1}) by \code{p2} (i.e. \code{m} stacked matrices each of \code{p1} by \code{p2})
#' or with dimensions \code{m} by (\code{p1} * \code{p2}) (i.e. regressor matrix with \code{m} regressor rows, each with column
#' length \code{p1 * p2}), respectively.
#' @param ridge_adj A ridge penalty tuning parameter (usually set to small value, default \code{1e-6}, to stabilize estimation).
#' @param wts the weight vector of \code{length(y)}. Default is 1.
#' @param int set to TRUE or FALSE to include intercept term in linear predictor (default \code{TRUE}).
#' @param M_pred (e.g. stacked (\code{q} * \code{p1}) by \code{p2}
#' signal inputs  or (unfolded) \code{q} by (\code{p1} * \code{p2}) signal
#' inputs for \code{q} new predictions.
#' @param  y_predicted a vector of responses from a cv data set (assoc. with \code{M_pred}), when \code{family = "gaussian"}.
#' @param family the response distribution, e.g.
#' \code{"gaussian", "binomial", "poisson", "Gamma"} distribution. Quotes are needed. Default is "gaussian".
#' @param m_binomial a vector of binomial trials having \code{length(y)}. Default is 1 vector for \code{family = "binomial"}, NULL otherwise.
#' @param r_gamma a vector of gamma shape parameters. Default is 1 vector for for \code{family = "Gamma"}, NULL otherwise.
#' @param link the link function, one of \code{"identity"}, \code{"log"}, \code{"sqrt"},
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, \code{"reciprocal"};
#' quotes are needed (default \code{"identity"}).
#' @param se_pred a scalar, e.g. \code{se = 2} (default) to produce twice se surfaces,
#' set \code{se} > 0. Used for CIs at \code{XYpred} locations.
#'
#' @return
#' \item{pcoef}{a vector of length \code{(Pars[1,3]+Pars[1,4])*(Pars[2,3]+Pars[2,4])}
#' of (unfolded) estimated P-spline coefficients for tensor surface.}
#' \item{summary_predicted}{inverse link prediction vectors, and standard error surfaces.}
#' \item{dev}{deviance of fit.}
#' \item{eff_df}{the approximate effective dimension of fit.}
#' \item{aic}{AIC.}
#' \item{df_resid}{approximate df residual.}
#' \item{cv}{leave-one-out standard error prediction, when \code{family = "gaussian"}.}
#' \item{cv_predicted}{standard error prediction for \code{y_predict}, when \code{family = "gaussian"}.}
#' \item{avediff_pred}{mean absolute difference prediction, when \code{family = 'gaussian'}.}
#' \item{Pars}{design and tuning parameters (see above arguments).}
#' \item{Dispersion_parm}{estimate of dispersion, \code{dev/df_resid}.}
#' \item{summary_predicted}{inverse link prediction vectors at \code{M_pred}, and standard error bands.}
#' \item{eta_predicted}{estimated linear predictor of \code{length(y)}.}
#' \item{press_mu}{leave-one-out prediction of mean, when \code{family = "gaussian"}.}
#' \item{bin_percent_correct}{percent correct classification based on 0.5 cut-off,
#' when \code{family = "binomial"}, NULL otherwise.}
#' \item{B}{Tensor basis (\code{p1} x \code{p2}) by (\code{n1} x \code{n2}) for 2D signal regression.}
#' \item{Q}{Effective regressors (\code{m} by \code{n1} * \code{n2}) for 2D signal regression.}
#' \item{Ahat}{smooth P-spline coefficient vector of length \code{p1} x \code{p2},
#' constructed by \code{B} \%*\% \code{pcoef}.}
#' \item{M}{the signal/image regressors.}
#' \item{y}{the response vector.}
#' \item{M1index}{index of length \code{p1} for rows of regressor matrix.}
#' \item{M2index}{index of length \code{p2} for columns of regressor matrix.}
#' \item{M_type}{"stacked" or "unfolded".}
#' \item{w}{GLM weight vector of length \code{m}.}
#' \item{h}{"hat" diagonals.}
#' \item{ridge_adj}{additional ridge tuning parameter to stabilize estimation.}



#' @author Paul Eilers and Brian Marx
#' @references Marx, B.D. and Eilers, P.H.C. (2005).
#' Multidimensional penalized signal regression, \emph{Technometrics}, 47: 13-22.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
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
#' # Plotting coefficient image
#'  plot(fit)
#' @export


"ps2DSignal" <-
  function(y, M, p1, p2, M_type = "stacked", M1_index = c(1:p1), M2_index = c(1:p2),
             Pars = rbind(c(1, p1, 10, 3, 1, 2), c(1, p2, 10, 3, 1, 2)), ridge_adj = 1e-6, M_pred = M, y_predicted = NULL,
             family = "gaussian", link = "default", m_binomial = 1 + 0*y, wts =
               1 + 0*y, r_gamma = 1 + 0*y, int = TRUE, se_pred = 2) {

    # Dimension, weight, parameter defaults
    m <- length(y)

    # 	Check to see if any argument is improperly defined
    parms <- pspline2d_checker(
      family, link, Pars[1, 4], Pars[2, 4],
      Pars[1, 6], Pars[2, 6], Pars[1, 3],
      Pars[2, 3], Pars[1, 5], Pars[2, 5],
      ridge_adj, wts
    )
    family <- parms$family
    link <- parms$link
    ridge_adj <- parms$ridge_adj
    wts <- parms$wts
    Pars[1, 3:6] <- c(
      parms$nseg1, parms$bdeg1, parms$lambda1,
      parms$pord1
    )
    Pars[2, 3:6] <- c(
      parms$nseg2, parms$bdeg2, parms$lambda2,
      parms$pord2
    )

    M <- X <- as.matrix(M)
    if (M_type == "stacked") {
      p1_ <- nrow(M) / m
      p2_ <- ncol(M)
      if (p1 != p1_ | p2 != p2_) {
        warning(paste(
          "recheck input p1 or p2: at least one dimension does not match"
        ))
      }
      x <- as.vector(t(M))
      X <- matrix(x, m, p1 * p2, byrow = TRUE)
    }

    # Prepare bases for estimation
    oM1 <- outer(rep(1, p2), M1_index)
    B1 <- bbase(as.vector(oM1), Pars[1, 1], Pars[1, 2], Pars[1, 3], Pars[1, 4])
    oM2 <- outer(M2_index, rep(1, p1))
    B2 <- bbase(as.vector(oM2), Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
    n1 <- ncol(B1)
    n2 <- ncol(B2)

    # Compute tensor products for estimated alpha surface
    B1_ <- kronecker(B1, t(rep(1, n2)))
    B2_ <- kronecker(t(rep(1, n1)), B2)
    B_ <- B1_ * B2_


    # Construct penalty matrices
    d1 <- Pars[1, 6]
    D1 <- diag(n1)
    D1 <- diff(D1, diff = d1)
    lambda1 <- Pars[1, 5]
    P1 <- sqrt(lambda1) * kronecker(D1, diag(n2))
    d2 <- Pars[2, 6]
    D2 <- diag(n2)
    D2 <- diff(D2, diff = d2)

    lambda2 <- Pars[2, 5]
    P2 <- sqrt(lambda2) * kronecker(diag(n1), D2)
    Pen <- rbind(P1, P2)

    # Additional ridge penalty construction to stablize estimation
    p_ridge <- NULL
    if (ridge_adj > 0) {
      nix_ridge <- rep(0, n1 * n2)
      p_ridge <- sqrt(ridge_adj) * diag(n1 * n2)
    }

    # Data augmentation and regression
    z1 <- rep(0, n2 * (n1 - d1))
    z2 <- rep(0, n1 * (n2 - d2))
    Q <- X %*% B_
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

    # Percent correct classification
    bin_percent_correct <- bin_0percent <- bin_1percent <- NULL
    if (family == "binomial") {
      count1 <- count2 <- pcount1 <- pcount2 <- 0
      p_hat <- mu / m_binomial
      for (ii in 1:m) {
        if (p_hat[ii] > 0.5) {
          count1 <- y[ii]
          count1 <- pcount1 + count1
          pcount1 <- count1
        }
        if (p_hat[ii] <= 0.5) {
          count2 <- m_binomial[ii] - y[ii]
          count2 <- pcount2 + count2
          pcount2 <- count2
        }
      }
      bin_percent_correct <- (count1 + count2) / sum(m_binomial)
      bin_1percent <- count1 / sum(m_binomial[p_hat > 0.5])
      bin_0percent <- count2 / sum(m_binomial[p_hat <= 0.5])
    }

    # GLM weights
    w <- ps_fit$w

    # Hat diagional useful for CV and ED
    qr_ <- ps_fit$f$qr
    h <- hat(qr_, intercept = FALSE)[1:m]
    eff_dim <- sum(h)

    # Deviance, by family
    dev_ = dev_calc(family, y, mu, m_binomial, r_gamma)
    dev = dev_$dev
    dispersion_parm = dev_$dispersion_parm

    # CV calculation, LOOCV
    cv <- press_mu <- press_e <- var_c <- NULL
    if (family == "gaussian") {
      dev <- sum(ps_fit$f$residuals[1:m]^2)
      dispersion_parm <- dev / (m - eff_dim)
      press_e <- ps_fit$f$residuals[1:m] / (1 - h)
      cv <- sqrt(sum((press_e)^2) / (m))
      press_mu <- y - press_e
    }

    # AIC
    aic <- dev + 2 * eff_dim
    w_aug <- c(w, (c(z1, z2) + 1))

    # Smooth coefficient surface
    if (int) {
      A_hat <- B_ %*% pcoef[2:(n_col + 1)]
      yint <- ps_fit$coef[1]
    }
    if (!int) {
      yint <- NULL
      A_hat <- B_ %*% pcoef
         }
    A_hatm <- matrix(A_hat, p1, p2, byrow = TRUE)

    # External prediction
    summary_predicted <- NULL
    cv_predicted <- eta_predicted <- avediff_pred <- NULL
    if (nrow(M_pred)>1) {
      q <- nrow(M_pred)
      X_p <- as.matrix(M_pred)
      if (M_type == "stacked") {
        M_pred <- as.matrix(M_pred)
        q <- nrow(M_pred) / p1
        x <- as.vector(t(M_pred))
        X_p <- matrix(x, q, p1 * p2, byrow = TRUE)
      }
      if (!int) {
        eta_predicted <- X_p %*% as.vector(A_hat)
        XB_ <- X_p %*% B_
      }
      if (int) {
        eta_predicted <- X_p %*% A_hat + yint
        XB_ <- cbind(rep(1, q), (X_p %*% B_))
      }

      # Covariances of coefficients
      R <- qr.R(qr_)
      # C2 = chol2inv(R)

      # Variances of fitted values
      L <- forwardsolve(t(R), t(XB_))
      var_pred <- colSums(L * L)
      stdev_pred <- as.vector(sqrt(var_pred))
      stdev_pred <- sqrt(dispersion_parm) * stdev_pred
      pivot <- as.vector(se_pred * stdev_pred)
      upper <- eta_predicted + pivot
      lower <- eta_predicted - pivot
      summary_predicted <- cbind(lower, eta_predicted, upper)
      if (!missing(y_predicted)) {
        if (family == "gaussian") {
          cv_predicted <- sqrt(sum((y_predicted -
            eta_predicted)^2) / (length(y_predicted)))
          avediff_pred <- (sum(y_predicted -
            eta_predicted)) / length(y_predicted)
        }
      }
      bin_percent_correct <- bin_0percent <- bin_1percent <- NULL
      if (link == "logit") {
        count1 <- count2 <- pcount1 <- pcount2 <- 0
        p_hat <- exp(eta_predicted) / (1 + exp(eta_predicted))
        if (!missing(y_predicted)) {
          y_predicted <- as.vector(y_predicted)
          for (ii in 1:length(eta_predicted)) {
            if (p_hat[ii] > 0.5) {
              count1 <- y_predicted[ii]
              count1 <- pcount1 + count1
              pcount1 <- count1
            }
            if (p_hat[ii] <= 0.5) {
              count2 <- 1 - y_predicted[ii]
              count2 <- pcount2 + count2
              pcount2 <- count2
            }
          }
          bin_percent_correct <- (count1 + count2) / length(
            y_predicted
          )
          bin_1percent <- count1 / sum(y_predicted)
          bin_0percent <- count2 / (length(y_predicted) -
            sum(y_predicted))
        }
      }

      summary_predicted <- inverse_link(summary_predicted, link)
      if (link == "reciprocal") {
        summary_predicted <- summary_predicted[, 3:1]
      }
      summary_predicted <- as.matrix(summary_predicted)
      dimnames(summary_predicted) <- list(NULL, c(
        "-se_pred * std_Lower",
        "Predicted", "+se_pred * std_Upper"
      ))
    }

    P <- list(
      pcoef = pcoef, Pars = Pars, yint = yint, int = int, family = family,
      link = link, dev = dev, aic = aic, bin_percent_correct =
      bin_percent_correct, bin_0 = bin_0percent, bin_1 = bin_1percent,
      df_resid = m - eff_dim, dispersion_parm = dispersion_parm, mu =
      mu, press_mu = press_mu, summary_predicted = summary_predicted,
      eta_predicted = eta_predicted, avediff_pred = avediff_pred,
      ridge_adj = ridge_adj, cv = cv, cv_predicted = cv_predicted,
      eff_dim = eff_dim, Q = Q, B = B_, h = h, Ahat = A_hat,
      M1index = M1_index, M2index = M2_index, y = y, M = M, M_type = M_type
    )
    class(P) <- "ps2dsignal"
    return(P)
  }
