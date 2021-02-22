#' Two-dimensional smoothing of scattered normal or non-normal (GLM)
#' responses using tensor product P-splines.
#'
#' @description \code{ps2DGLM} is used to smooth scattered
#' normal or non-normal responses, with aniosotripic
#' penalization of tensor product P-splines.
#'
#'
#' @details Support functions needed: \code{pspline_fitter}, \code{bbase}, and \code{pspline_2dchecker}.
#' @seealso ps2DNormal
#' @import stats
#'
#' @param Data a matrix of 3 columns \code{x, y, z} of equal length;
#' the response is \code{z}.
#' @param Pars a matrix of 2 rows, where the first and second row
#' sets the P-spline paramters for \code{x} and \code{y}, respectively.
#' Each row consists of: \code{min max nseg bdeg lambda pord}.
#' The \code{min} and \code{max} set the ranges, \code{nseg} (default 10)
#' is the number of evenly spaced segments between \code{min} and \code{max},
#' \code{bdeg} is the degree of the basis (default 3 for cubic),
#' \code{lambda} is the (positive) tuning parameter for the penalty (default 1),
#' \code{pord} is the number for the order of the difference penalty (default 2).
#' @param XYpred a matrix with two columns \code{(x, y)} that give the coordinates
#' of (future) prediction; the default is the data locations.
#' @param family \code{"gaussian", "binomial", "poisson", "Gamma"} (quotes needed). Default is "gaussian".
#' @param link the link function, one of \code{"identity"}, \code{"log"}, \code{"sqrt"},
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, \code{"reciprocal"};
#' quotes are needed (default \code{"identity"}).
#' @param wts  non-negative weights, which can be zero (default ones).
#' @param m_binomial vector of binomial trials, default is vector of ones with \code{family = "binomial"}, NULL otherwise.
#' @param r_gamma gamma scale parameter, default is vector ones with \code{family = "Gamma"}, NULL otherwise.
#' @param ridge_adj a ridge penalty tuning parameter, usually set to small value, e.g. \code{1e-8} to stabilize estimation (default 0).
#' @param se_pred a scalar, default \code{se_pred = 2} to produce se surfaces,
#' set \code{se_pred} > 0. Used for CIs for \code{XYpred} locations.
#' @param z_predicted a vector of responses associated with \code{XYpred}, useful for external validation with \code{family = "gaussian"}.
#'
#' @return
#' \item{pcoef}{a vector of length \code{(Pars[1,3]+Pars[1,4])*(Pars[2,3]+Pars[2,4])}
#' of (unfolded) estimated P-spline coefficients.}
#' \item{mu}{a vector of \code{length(z)} of smooth estimated means (at the \code{x,y} locations).}
#' \item{dev}{the deviance of fit.}
#' \item{eff_df}{the approximate effective dimension of fit.}
#' \item{aic}{AIC.}
#' \item{df_resid}{approximate df residual.}
#' \item{cv}{leave-one-out standard error prediction, when \code{family = 'gaussian'}.}
#' \item{cv_predicted}{standard error prediction for \code{y_predict}, when \code{family = 'gaussian'}.}
#' \item{avediff_pred}{mean absolute difference prediction, when \code{family = 'gaussian'}.}
#' \item{Pars}{the design and tuning parameters (see arguments above).}
#' \item{dispersion_parm}{estimate of dispersion, \code{dev/df_resid}.}
#' \item{summary_predicted}{inverse link prediction vectors, and \code{se_pred} bands.}
#' \item{eta_predicted}{estimated linear predictor of \code{length(z)}.}
#' \item{press_mu}{leave-one-out prediction of mean, when \code{family = 'gaussian'}.}
#' \item{bin_percent_correct}{percent correct classification based on 0.5 cut-off (when \code{family = "binomial"}).}
#' \item{Data}{a matrix of 3 columns \code{x, y, z} of equal length;
#' the response is \code{z}.}
#' \item{Q}{the tensor product B-spline basis.}
#' \item{qr}{the Q-R of the model.}
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
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
#' plot(fit, xlab = "Start", ylab = "Age")
#' #title(main = "Probability of Kyphosis")
#' @export

"ps2DGLM" <-
  function(Data, Pars = rbind(c(min(Data[, 1]), max(Data[, 1]), 10, 3, 1, 2),
                              c(min(Data[, 2]), max(Data[, 2]), 10, 3, 1, 2)),
                 ridge_adj = 0, XYpred = Data[, 1:2], z_predicted = NULL,
                 se_pred = 2, family = "gaussian", link = "default",
                 m_binomial = rep(1, nrow(Data)), wts = rep(1, nrow(Data)), r_gamma = rep(1, nrow(Data))) {

    # Prepare bases for estimation
    z <- Data[, 3]
    x <- Data[, 1]
    y <- Data[, 2]
    m <- length(z)

    # Check to see if any argument is improperly defined
    parms <- pspline2d_checker(
      family, link, Pars[1, 4], Pars[2, 4], Pars[1, 6],
      Pars[2, 6], Pars[1, 3], Pars[2, 3], Pars[1, 5], Pars[2, 5],
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
    B1 <- bbase(as.vector(x), Pars[1, 1], Pars[1, 2], Pars[1, 3], Pars[1, 4])
    B2 <- bbase(as.vector(y), Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
    n1 <- ncol(B1)
    n2 <- ncol(B2)

    # Compute tensor products for estimated alpha surface
    B1_ <- kronecker(B1, t(rep(1, n2)))
    B2_ <- kronecker(t(rep(1, n1)), B2)
    Q <- B1_ * B2_

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

    # Ridge penalty
    p_ridge <- 0*diag(n1*n2)
    if (ridge_adj > 0) {
      nix_ridge <- rep(0, n1 * n2)
      p_ridge <- sqrt(ridge_adj) * diag(n1 * n2)
    }

    # Data augmentation and regression
    z1 <- rep(0, n2 * (n1 - d1))
    z2 <- rep(0, n1 * (n2 - d2))
    n_col <- ncol(Q)

    # Apply the penalized GLM fitter on for y on Q
    ps_fit <- pspline_fitter(y = z, B = Q, family, link,
      P = Pen, P_ridge = p_ridge, wts, m_binomial, r_gamma)

    mu <- ps_fit$mu
    pcoef <- ps_fit$coef

    # Percent correct classification for binomial
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

    # Hat diagaonal for ED
    qr <- ps_fit$f$qr
    h <- hat(qr, intercept = FALSE)[1:m]
    eff_dim <- sum(h)

    # Deviance, by family
    dev_ = dev_calc(family, y = z, mu, m_binomial, r_gamma)
    dev = dev_$dev
    dispersion_parm = dev_$dispersion_parm

    # Leave-one-out CV for Gaussian	response
    cv <- press_mu <- press_e <- var_c <- NULL

    if (family == "gaussian") {
      dev <- sum(ps_fit$f$residuals[1:m]^2)
      dispersion_parm <- dev / (m - eff_dim)
      press_e <- ps_fit$f$residuals[1:m] / (1 - h)
      cv <- sqrt(sum((press_e)^2) / (m))
      press_mu <- z - press_e
    }
    aic <- dev + 2 * eff_dim

    # Weight augmentation
    w_aug <- c(w, (c(z1, z2) + 1))

    # Compute XY new predictions, and CIs
    summary.predicted = NULL
    cv_predicted <- eta_predicted <- avediff_pred <- NULL
    if (nrow(XYpred)>1) {
      Bxp <- bbase(XYpred[, 1], Pars[1, 1], Pars[1, 2], Pars[1, 3], Pars[1, 4])
      Byp <- bbase(XYpred[, 2], Pars[2, 1], Pars[2, 2], Pars[2, 3], Pars[2, 4])
      B1p <- kronecker(Bxp, t(rep(1, n2)))
      B2p <- kronecker(t(rep(1, n1)), Byp)
      Bp <- B1p * B2p
      eta_predicted <- Bp %*% pcoef

      # Covariances of coefficients
      R <- qr.R(qr)
      # Variances of fitted values
      L <- forwardsolve(t(R), t(Bp))
      var_pred <- colSums(L * L)
      stdev_pred <- as.vector(sqrt(var_pred))
      stdev_pred <- sqrt(dispersion_parm) * stdev_pred
      pivot <- as.vector(se_pred * stdev_pred)
      upper <- eta_predicted + pivot
      lower <- eta_predicted - pivot
      ieta_predicted <- inverse_link(eta_predicted, link)
      ilower <- inverse_link(lower, link)
      iupper <- inverse_link(upper, link)
      if (link == "reciprocal") {
        ilowup <- cbind(ilower, iupper)
        iupper <- ilowup[, 1]
        ilower <- ilowup[, 2]
      }
      summary_predicted <- as.matrix(cbind(
        ilower, ieta_predicted,
        iupper
      ))
      if (!missing(z_predicted)) {
        if (family == "gaussian") {
          cv_predicted <- sqrt(sum((z_predicted -
            eta_predicted)^2) / (length(z_predicted)))
          avediff_pred <- (sum(z_predicted -
            eta_predicted)) / length(z_predicted)
        }
      }
      bin_percent_correct <- bin_0percent <- bin_1percent <- NULL
      if (link == "logit") {
        count1 <- count2 <- pcount1 <- pcount2 <- 0
        p_hat <- exp(eta_predicted) / (1 + exp(eta_predicted))
        if (!missing(z_predicted)) {
          z_predicted <- as.vector(z_predicted)
          for (ii in 1:length(eta_predicted)) {
            if (p_hat[ii] > 0.5) {
              count1 <- z_predicted[ii]
              count1 <- pcount1 + count1
              pcount1 <- count1
            }
            if (p_hat[ii] <= 0.5) {
              count2 <- 1 - z_predicted[ii]
              count2 <- pcount2 + count2
              pcount2 <- count2
            }
          }
          bin_percent_correct <- (count1 + count2) / length(
            z_predicted
          )
          bin_1percent <- count1 / sum(z_predicted)
          bin_0percent <- count2 / (length(z_predicted) -
            sum(z_predicted))
        }
      }
      dimnames(summary_predicted) <- list(NULL, c(
        "-2std_Lower",
        "Predicted", "+2std_Upper"
      ))
    }
    P <- list(
      pcoef = pcoef, Pars = Pars, family = family, link = link, dev
      = dev, aic = aic, bin_percent_correct = bin_percent_correct,
      bin_0 = bin_0percent, bin_1 = bin_1percent, df_resid = m -
        eff_dim, dispersion_parm = dispersion_parm, mu = mu, press_mu =
        press_mu, summary_predicted = summary_predicted, eta_predicted
      = eta_predicted, avediff_pred = avediff_pred, ridge_adj =
        ridge_adj, cv = cv, cv_predicted = cv_predicted, eff_dim =
        eff_dim, Data = Data, Q = Q, qr = qr)
    class(P) <- "ps2dglm"
    return(P)
  }
