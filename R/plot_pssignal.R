#' Plotting function for \code{psSignal}
#'
#' @description Plotting function for signal regression P-spline smooth coefficients (using \code{psSignal} with \code{class pssignal}), with or
#' without standard error bands.
#'
#' @import graphics
#'
#' @param x the P-spline x, usually from \code{psSignal}.
#' @param ... other parameters.
#' @param se a scalar, e.g. \code{se = 2} to produce twice se bands, set \code{se} > 0 (or set \code{se = 0} to supress).
#' @param xlab label for the x-axis, e.g. "my x" (quotes required).
#' @param ylab label for the y-axis, e.g. "my y" (quotes required).
#' @param col color.
#' @param lty line type for plotting e.g. \code{lty = 2}.

#' @return
#' \item{Plot}{a plot of the smooth P-spline signal coefficent vector, with or without standard error bands.}
#'
#' @author  Paul Eilers and Brian Marx
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
#' X=X[sel, ]
#' iindex=iindex[sel]
#' dX=diff(X)
#' diindex=iindex[-1]
#' y=as.vector(labc[1,1:40])
#' oout = 23
#' dX=t(dX[,-oout])
#' y=y[-oout]
#' fit2 = psSignal(y, dX, diindex, nseg = 25,lambda = 0.0001)
#' plot(fit2, se = 2, xlab = 'Coefficient Index', ylab= "ps Smooth Coeff")
#' title(main='25 B-spline segments with tuning=0.0001')
#' names(fit2)
#'
#' @export
plot.pssignal = function(x,..., se = 2, xlab="", ylab="", col='black', lty=1){
  pssig = x
  beta=as.vector(pssig$beta)
  x_index=pssig$x_index
  if(se == 0){
      plot(x_index, beta, type = "l", lty = lty, col = col,
      xlab = xlab, ylab = ylab)}

    if (se > 0) {
      pivot=as.vector(se*pssig$stdev_beta)
      ubeta=beta+pivot
      lbeta=beta-pivot
      plot(x_index, beta, type = "l", lty = lty, col = col,
           xlab = xlab, ylab = ylab,
           ylim=c(min(lbeta),max(ubeta)))
      lines(x_index, ubeta, col = 2, type = "l", lty = 2)
      lines(x_index, lbeta, col = 2, type = "l", lty = 2)
    }
     if (se < 0) {
    warning(paste("se should be nonnegative"))
  }
}
