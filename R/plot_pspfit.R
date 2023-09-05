#' Plotting function for \code{psNormal}, \code{psPoisson}, \code{psBinomial}
#'
#' @description Plotting function for P-spline smooth with normal, Poisson, or binomial responses
#' (\code{class pspfit}), with or without standard error bands.
#'
#' @import graphics
#'
#' @param x the P-spline object, usually from psNormal, psPoisson, psBinomial.
#' @param ... other parameters.
#' @param se a scalar, e.g. \code{se = 2} to produce twice se bands, set \code{se} > 0 (or set \code{se}=0 to supress).
#' @param xlab label for the x-axis.
#' @param ylab label for the y-axis.
#' @param col color for points.
#' @param pch point character.

#' @return
#' \item{Plot}{a plot of the mean (inverse link) smoothed normal, Poisson, or binomial responses, with or without se bands.}
#'
#' @author  Paul Eilers and Brian Marx
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references  Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#' @import graphics
#' @examples
#' library(JOPS)
#' #Extract data
#' library(MASS)
#' # Get the data
#' data(mcycle)
#' x = mcycle$times
#' y = mcycle$accel
#' fit1 = psNormal(x, y, nseg = 20, bdeg = 3, pord = 2, lambda = .8)
#' plot(fit1, se = 2, xlab = "time (ms)", ylab = "accel")
#'
#' @examples
#' library(JOPS)
#' library(boot)
#' # Extract the data
#' Count = hist(boot::coal$date, breaks=c(1851:1963), plot = FALSE)$counts
#' Year = c(1851:1962)
#' xl = min(Year)
#' xr = max(Year)
#'
#' # Poisson smoothing
#' nseg = 20
#' bdeg = 3
#' fit1=psPoisson(Year, Count, xl, xr, nseg, bdeg, pord = 2,
#' lambda = 1)
#' names(fit1)
#' plot(fit1, xlab = "Year", ylab = "Count", se = 2)
#'
#' @examples
#' library(JOPS)
#' #Extract data
#' library(rpart)
#' Kyphosis = kyphosis$Kyphosis
#' Age  =kyphosis$Age
#' y = 1 * (Kyphosis == "present")  # make y 0/1
#' # Binomial smoothing
#' fit1 = psBinomial(Age, y, xl = min(Age), xr = max(Age), nseg = 20,
#'                  bdeg = 3, pord = 2, lambda = 1)
#' names(fit1)
#' plot(fit1, xlab = "Age", ylab = '0/1', se = 2)
#'
#' @export

plot.pspfit = function(x,..., se = 2, xlab = "", ylab = "", col='black', pch=1){

# Retrieve gridded x and linear predictor
ps = x
xgrid = ps$xgrid
ygrid = ps$ygrid
se_eta = ps$se_eta


# Code block for psNormal
  if (ps$family == 'gaussian') {
    plot(ps$x, ps$y, main = "", xlab = xlab,
         ylab = ylab, pch = pch, col = col)
    lines(xgrid, ygrid, col = 'blue')

    # Error bands (Bayesian estimate)
    if (se > 0) {
      lgrid = ygrid - se * se_eta
      ugrid = ygrid + se * se_eta
      lines(xgrid, lgrid, lty = 2, col = 'red')
      lines(xgrid, ugrid, lty = 2, col = 'red')
    }
  }
# Code block for psPoisson
  if(ps$family=='poisson'){
    # Compute curve on grid
    mugrid = exp(ygrid)
    # Plot data and fit
    plot(ps$x, ps$y, xlab = xlab, ylab = ylab,
         col = col, pch = pch)
    lines(xgrid, mugrid, col = "blue")

    # SE bands on a grid
    if(se>0){
      ugrid = exp(ygrid + se*se_eta)
      lgrid = exp(ygrid - se*se_eta)
      lines(xgrid, ugrid, col='red')
      lines(xgrid, lgrid, col='red')
    }
  }
# Code block for psBinomial
  if(ps$family == "binomial"){
    # Compute curve on grid
    # Plot data and fit
    plot(ps$x, ps$y, xlab = xlab, ylab = ylab,
         col = col, pch = pch)
    lines(xgrid, ps$pgrid, col = "blue")
  # SE bands on a grid
    if(se>0){
      ugrid = exp(ygrid + se*se_eta)/(1 + exp(ygrid +  se*se_eta) )
      lgrid = exp(ygrid - se*se_eta)/(1+exp(ygrid - se*se_eta) )
      lines(xgrid, lgrid, col='red')
      lines(xgrid, ugrid, col='red')
    }
  }
    if (se < 0) {
    warning(paste("se should be nonnegative"))
  }
}
