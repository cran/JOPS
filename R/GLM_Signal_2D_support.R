#' P-spline fitting algorithm for the GLM.
#'
#' @description \code{pspline_fitter} appies the method of scoring
#' to a variety of response distributions and link functions within
#' for P-spline fitting within the GLM framework.
#'
#' @param y the glm response vector of length \code{m}.
#' @param B The effective P-spline regressors, e.g. \code{B} for B-splines, \code{Q=X \%*\% B} for PSR.
#' @param family the response distribution, e.g. \code{"gaussian", "binomial", "poisson", "Gamma"} distribution; quotes are needed
#' (default \code{family = "gaussian"}.)
#' @param link the link function, one of \code{"identity"}, \code{"log"}, \code{"sqrt"},
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, \code{"reciprocal"};
#' quotes are needed (default \code{link = "identity"}).
#' @param P P-spline ("half") penalty matrix for data augmentation, such that \code{P'P = lambda D'D}.
#' @param P_ridge ridge ("half") penalty for data augmentation, usually \code{sqrt(lambda_r)*I} (default 0).
#' @param wts the weight vector of \code{length(y)}, separate from GLM weights.
#' @param m_binomial a vector of binomial trials having \code{length(y)}, when \code{family = "binomial"}.
#' Default is 1 vector.
#' @param r_gamma a vector of gamma shape parameters, when \code{family = "Gamma"}. Default is 1 vector.

#' @return
#' \item{coef}{the estimated P-spline coefficient regressor, using the effective regressors.}
#' \item{w}{\code{wts*w}, GLM weight vector times input weights of length m.}
#' \item{f}{the \code{lsfit} object using data augmentation to get P-spline coefficient estimates.}
#' \item{eta}{the linear predictor from \code{f}.}

#' @export
"pspline_fitter" <-
  function(y, B, family="gaussian", link="identity", P, P_ridge = 0*diag(ncol(B)),
     wts = 0*y + 1, m_binomial = 0*y + 1, r_gamma = 0*y + 1)
  {
    coef_est <- rep(1, ncol(B))
    sumpr = sum(P_ridge)
    if(sumpr > 0){
      nix_ridge = rep(0, nrow(P_ridge))
    }
    nix = rep(0, nrow(P))
    if (family == "binomial") {
      mu <- (y + 0.5 * m_binomial) / 2
    }
    if (family == "Gamma" || family == "poisson") {
      mu <- (y + 3)
    }
    if (family == "gaussian") {
      mu <- rep(mean(y), length(y))
    }
    if (link == "identity") {
      eta <- mu
    }
    if (link == "log") {
      eta <- log(mu)
    }
    if (link == "sqrt") {
      eta <- sqrt(mu)
    }
    if (link == "logit") {
      eta <- log(mu / (m_binomial - mu))
    }
    if (link == "reciprocal") {
      eta <- 1 / mu
    }
    if (link == "probit") {
      eta <- qnorm(mu / m_binomial)
    }
    if (link == "cloglog") {
      eta <- log(-log(1 - mu / m_binomial))
    }
    if (link == "loglog") {
      eta <- -log(-log(mu / m_binomial))
    }
    for (ii in 1:100) {
      if (ii > 100) {
        break
      }
      if (link == "identity") {
        mu <- eta
        h_prime <- 1
      }
      if (link == "log") {
        mu <- exp(eta)
        h_prime <- mu
      }
      if (link == "sqrt") {
        mu <- eta^2
        h_prime <- 2 * eta
      }
      if (link == "logit") {
        mu <- m_binomial / (1 + exp(-eta))
        h_prime <- mu * (1 - mu / m_binomial)
      }
      if (link == "reciprocal") {
        mu <- 1 / eta
        h_prime <- -(mu^2)
      }
      if (link == "probit") {
        mu <- m_binomial * pnorm(eta)
        h_prime <- m_binomial * dnorm(eta)
      }
      if (link == "cloglog") {
        mu <- m_binomial * (1 - exp(-exp(eta)))
        h_prime <- (m_binomial) * exp(eta) * exp(-exp(eta))
      }
      if (link == "loglog") {
        mu <- m_binomial * exp(-exp(-eta))
        h_prime <- m_binomial * exp(-eta) * exp(-exp(-eta))
      }
      if (family == "gaussian") {
        w <- rep(1, length(y))
      }
      if (family == "poisson") {
        w <- h_prime^2 / mu
      }
      if (family == "binomial") {
        w <- h_prime^2 / (mu * (1 - mu / m_binomial))
      }
      if (family == "Gamma") {
        w <- (r_gamma * h_prime^2) / mu^2
      }
      u <- (y - mu) / h_prime + eta
      if (sumpr > 0) {
        f <- lsfit(rbind(B, P, P_ridge), c(u, nix, nix_ridge),
          wt = c(wts, nix + 1, nix_ridge + 1) * c(w, (nix +
            1), (nix_ridge + 1)), intercept = FALSE)
      }
      if (sumpr == 0) {
        f <- lsfit(rbind(B, P), c(u, nix), wt = c(wts, nix + 1) *
          c(w, (nix + 1)), intercept = FALSE)
      }
      coef_old <- coef_est
      coef_est <- as.vector(f$coef)
      eta <- B %*% coef_est
      if (link == "identity") {
        mu <- eta}
      d_coef <- max(abs((coef_est - coef_old) / coef_old))
      if(family == 'gaussian' && link == 'identity'){
              d_coef = 0
              }
      print(c(ii, d_coef))
      if (d_coef < 1e-006) {
        break
      }
    }
    if (ii > 99) {
      warning(paste("parameter estimates did NOT converge in 100 iterations"))
    }
    llist <- list(coef = coef_est, mu = mu, f = f, w = w * wts, eta = eta)
    return(llist)

    }

#' P-spline checking algorithm for the GLM.
#'
#' @description \code{pspline_checker} checks to see if all the inputs associated
#' for P-spines are properly defined.
#'
#' @param family the response distribution, e.g. \code{"gaussian", "binomial", "poisson", "Gamma"} distribution. Quotes are needed.
#' @param link the link function, one of \code{"identity"}, \code{"log"}, \code{"sqrt"},
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, \code{"reciprocal"};
#' @param bdeg the degree of B-splines.
#' @param pord the order of the penalty.
#' @param nseg the number of evenly-spaced B-spline segmements.
#' @param wts the weight vector, separate from GLM weights.
#' @param lambda the positive tuning parameter for the difference penalty.
#' @param ridge_adj the positive tuning parameter for the ridge penalty.
#'
#' @return
#' \item{list}{same as inputs, with warnings if required.}
#' @export
"pspline_checker" <-
  function(family, link, bdeg, pord, nseg, lambda, ridge_adj, wts) {
    if (link == "default" && family == "gaussian") {
      link <- "identity"
    }
    if (link == "default" && family == "poisson") {
      link <- "log"
    }
    if (link == "default" && family == "binomial") {
      link <- "logit"
    }
    if (link == "default" && family == "Gamma") {
      link <- "log"
    }
    if (family != "binomial" && family != "gaussian" && family != "poisson" &&
        family != "Gamma") {
      warning(paste("Improper FAMILY option. Choose: gaussian, poisson, binomial or Gamma"))
    }
    if ((family == "binomial") && (link != "logit" && link != "probit" &&
                                   link != "cloglog" && link != "loglog")) {
      warning(paste("Improper LINK option with family=binomial. Choose: logit, probit, loglog, cloglog"))
    }
    if ((family == "Gamma") && (link != "log" && link != "reciprocal" && link !=
                                "identity")) {
      warning(paste("Improper LINK option with family=Gamma. Choose: reciprocal, log, identity"))
    }
    if ((family == "poisson") && (link != "log" && link != "sqrt" && link !=
                                  "identity")) {
      warning(paste("Improper LINK option with family=poisson. Choose: log, sqrt, identity"))
    }
    if ((family == "gaussian") && (link != "identity")) {
      warning(paste("Improper LINK option with family=gaussian. Choose: identity"))
    }
    if (bdeg < 0) {
      bdeg <- 1
      warning(paste("bdeg must be non-neg integer: have used 1"))
    }
    if (pord < 0) {
      pord <- 0
      warning(paste("pord must be non-neg integer: have used 0"))
    }
    if (nseg < 2) {
      nseg <- 2
      warning(paste("nseg must be positive integer, > 1: have used 2"))
    }
    if (lambda < 0) {
      lambda <- 0
      warning(paste("lambda cannot be negative: have used 0"))
    }
    if (ridge_adj < 0) {
      ridge_adj <- 0
      warning(paste("ridge_adj cannot be negative: have used 0"))
    }
    if (min(wts) < 0) {
      warning(paste("At least one weight entry is negative"))
    }
    llist <- list(
      family = family, link = link, bdeg = bdeg, pord =
        pord, nseg = nseg, lambda = lambda, ridge_adj
      = ridge_adj, wts = wts
    )
    return(llist)
  }

#' P-spline 2D tensor product checking algorithm for the GLM.
#'
#' @description \code{pspline_2dchecker} checks to see if all the 2D tensor inputs associated
#' for P-spines are properly defined.

#'
#' @param family the response distribution, e.g. \code{"gaussian", "binomial", "poisson", "Gamma"} distribution. Quotes are needed.
#' @param link the link function, one of \code{"identity"}, \code{"log"}, \code{"sqrt"},
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, \code{"reciprocal"};
#' quotes are needed.
#' @param bdeg1 the degree of B-splines.
#' @param bdeg2 the degree of B-splines.
#' @param pord1 the order of the penalty.
#' @param pord2 the order of the penalty.
#' @param nseg1 the number of evenly spaced B-spline segmements.
#' @param nseg2 the number of evenly spaced B-spline segmements.
#' @param wts the weight vector, separate from GLM weights.
#' @param lambda1 the positive tuning parameter for the difference penalty.
#' @param lambda2 the positive tuning parameter for the difference penalty.
#' @param ridge_adj the positive tuning parameter for the ridge penalty.
#'
#' @return
#' \item{list}{same as inputs, with warnings if required.}

#' @export
"pspline2d_checker" <-
  function(family, link, bdeg1, bdeg2, pord1, pord2, nseg1,
           nseg2, lambda1, lambda2, ridge_adj, wts) {
    if (link == "default" && family == "gaussian") {
      link <- "identity"
    }
    if (link == "default" && family == "poisson") {
      link <- "log"
    }
    if (link == "default" && family == "binomial") {
      link <- "logit"
    }
    if (link == "default" && family == "Gamma") {
      link <- "log"
    }
    if (family != "binomial" && family != "gaussian" && family != "poisson" &&
        family != "Gamma") {
      warning(paste("Improper FAMILY option. Choose: gaussian, poisson, binomial or Gamma"))
    }
    if ((family == "binomial") && (link != "logit" && link != "probit" &&
                                   link != "cloglog" && link != "loglog")) {
      warning(paste("Improper LINK option with family=binomial. Choose: logit, probit, loglog, cloglog"))
    }
    if ((family == "Gamma") && (link != "log" && link != "reciprocal" && link !=
                                "identity")) {
      warning(paste("Improper LINK option with family=Gamma. Choose: reciprocal, log, identity"))
    }
    if ((family == "poisson") && (link != "log" && link != "sqrt" && link !=
                                  "identity")) {
      warning(paste("Improper LINK option with family=poisson. Choose: log, sqrt, identity"))
    }
    if ((family == "gaussian") && (link != "identity")) {
      warning(paste("Improper LINK option with family=gaussian. Choose: identity"))
    }
    if (bdeg1 < 0) {
      bdeg1 <- 1
      warning(paste("bdeg1 must be non-neg integer: have used 1"))
    }
    if (pord1 < 0) {
      pord1 <- 0
      warning(paste("pord1 must be non-neg integer: have used 0"))
    }
    if (nseg1 < 2) {
      nseg1 <- 2
      warning(paste("nseg1 must be positive integer, > 1: have used 2"))
    }
    if (lambda1 < 0) {
      lambda1 <- 0
      warning(paste("lambda1 cannot be negative: have used 0"))
    }
    if (bdeg2 < 0) {
      bdeg2 <- 1
      warning(paste("bdeg2 must be non-neg integer: have used 1"))
    }
    if (pord2 < 0) {
      pord2 <- 0
      warning(paste("pord2 must be non-neg integer: have used 0"))
    }
    if (nseg2 < 2) {
      nseg2 <- 2
      warning(paste("nseg2 must be positive integer, > 1: have used 2"))
    }
    if (lambda2 < 0) {
      lambda2 <- 0
      warning(paste("lambda2 cannot be negative: have used 0"))
    }
    if (ridge_adj < 0) {
      ridge_adj <- 0
      warning(paste("ridge_adj cannot be negative: have used 0"))
    }
    if (min(wts) < 0) {
      warning(paste("At least one weight entry is negative"))
    }
    llist <- list(
      family = family, link = link, bdeg1 = bdeg1, pord1
      = pord1, nseg1 = nseg1, lambda1 = lambda1,
      bdeg2 = bdeg2, pord2 = pord2, nseg2 =
        nseg2, lambda2 = lambda2, ridge_adj = ridge_adj, wts =
        wts
    )
    return(llist)
  }



#' Deviance calculation for GLM P-spline fitting.
#'
#' @description Calculates the deviance and returns the ML estimated dispersion parameter
#' for a variety of response distributions for P-spline fitting within the GLM framework.
#'
#' @param family the response distribution, e.g. \code{"gaussian", "binomial", "poisson", "Gamma"} distribution. Quotes are needed; default \code{"family = gaussian"}.
#' @param y the glm response vector of length \code{m}.
#' @param mu the P-spline estimated mean for the glm response vector of length \code{m}.
#' @param m_binomial a vector of binomial trials having \code{length(y)}, when \code{family = "binomial"}. Default is 1 vector.
#' @param r_gamma a vector of gamma shape parameters, when \code{family = "Gamma"}. Default is 1 vector.

#' @return A list with two fields:
#' \item{dev}{the estimated deviance.}
#' \item{dispersion_parm}{the ML estimated dispersion parameter.}

#' @export

dev_calc = function(family = "gaussian", y, mu, m_binomial =
                      0*y + 1, r_gamma = 0*y + 1 ) {
  if (family == "gaussian") {
    dev <- sum((y-mu)^2)
    dispersion_parm = dev/length(y)
  }
  e = 1e-9
  if (family == "binomial") {
  dev <- 2 * sum((y + e) * log((y + e) / mu) + (m_binomial - y + e) *
                   log((m_binomial - y + e) / (m_binomial - mu)))
  dispersion_parm <- 1
  }
if (family == "poisson") {
  dev <- 2 * sum(y * log(y + e) - y - y * log(mu) + mu)
  dispersion_parm <- 1
  }
if (family == "Gamma") {
  dev <- -2 * sum(r_gamma * (log((y + e) / mu) - ((y - mu) / mu)))
  ave_dev <- dev / length(y)
  dispersion_parm <- (ave_dev * (6 + ave_dev)) / (6 + 2 * ave_dev)
  }
llist = list(dev = dev, dispersion_parm = dispersion_parm)
return(llist)
}
