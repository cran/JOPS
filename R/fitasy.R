#' Fit asymmetry parameters in the expectile bundle model
#'
#' @description There are two functions for fitting the expectile bundle model, the present one for estimating asymmetry parameters (\code{fitasy}),
#' the other for estimating the amplitude function, \code{fitampl}. See the details below.
#'
#' @param y a response vector.
#' @param B a proper B-spline basis matrix, see \code{bbase()}.
#' @param b a vector of B-spline coefficients.
#' @param p a vector of asymmetries with values between  0 and 1.
#' @param c0 a vector.
#'
#' @return a vector of estimated asymmetry parameters .
#'
#' @details
#' The expectile bundle model determines a set of expectile curves for a point cloud with data vectors \code{x} and \code{y},
#' as \eqn{\psi_j{x_i} = a_j g(x_i)}. Here \eqn{a_j} is the asymmetry parameter corresponding to a given asymmetry \eqn{p_j}.
#' A vector of asymmetries with all \eqn{0 <p_j < 1} is specified by the user.
#'
#' The asymmetric least squares objective function is
#' \deqn{\sum_j \sum_i w_{ij}(y_i - \sum_j a_j g_j(x_i))^2.}
#' The function \eqn{g(\cdot)} is called the amplitude. The weights depend on the residuals:
#' \deqn{w_{ij} = p_j} if \eqn{y_i > a_jg(x_i)} and \eqn{w_{ij} = 1- p_j} otherwise.
#'
#' The amplitude function is a sum of B-splines with coefficients \code{alpha}. There is no direct solution, so \code{alpha}
#' and the asymmetry parameters \code{a} must be updated alternatingly. See the example.
#'
#' @note  This is a simplification of the model described in the reference. There is no explict term for the trend.
#'
#' @author Paul Eilers
#'
#' @references Schnabel, S.K. and Eilers, P.H.C. (2013) A location-scale model for non-crossing expectile curves. \emph{Stat} 2: 171–183.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @examples
#' # Get the data
#' data(bone_data)
#' x = bone_data$age
#' y = bone_data$spnbmd
#' m <- length(x)
#'
#' # Set asymmetry levels
#' p = c(0.005, 0.01, 0.02, 0.05, 0.2, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995)
#' np <- length(p)
#'
#' # Set P-spline parameters
#' x0 <- 5
#' x1 <- 30
#' ndx <- 20
#' bdeg <- 3
#' pord <- 2
#'
#' # Compute bases
#' B <- bbase(x, x0, x1, ndx, bdeg)
#' xg <- seq(from = min(x), to = max(x), length = 100)
#' Bg <- clone_base(B, xg)
#' n <- ncol(B)
#'
# Fit the model
#' lambda = 1
#' alpha <- rep(1,n)
#' a = p
#' for (it in 1:20){
#'   alpha <- fitampl(y, B, alpha, p, a, pord, lambda)
#'   alpha <- alpha / sqrt(mean(alpha ^ 2))
#'   anew <- fitasy(y, B, alpha, p, a)
#'   da = max(abs(a - anew))
#'   a = anew
#'   cat(it, da, '\n')
#'      if (da < 1e-6) break
#' }
#'
#' # Compute bundle on grid
#' ampl <- Bg %*% alpha
#' Z <- ampl %*% a
#'
#' # Plot data and bundle
#' plot(x, y, pch = 15, cex = 0.7, col = 'grey', xlab = 'Age', ylab = 'Density')
#' cols = colorspace::rainbow_hcl(np, start = 10, end = 350)
#' matlines(xg, Z, lty = 1, lwd = 2, col = cols)
#'
#'
#' @export

fitasy <- function(y, B, b, p, c0){
  a <- B %*% b
  ccr <- 0 * p
  for(j in 1:length(p)){
    w <- ifelse(y >= c0[j]*a, p[j], 1 - p[j])
    ccr[j] <- sum(w * a * y) / sum(w * a * a)
  }
  return(ccr)
}

