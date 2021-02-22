#' Prepare graphics layout for multiple panels
#'
#' @description Adapt margins and axes layout for multiple panels.
#'
#' @param rows number of rows.
#' @param cols number of columns.
#'
#' @return Prepare graphics layout for multiple panels
#'
#' @author Paul Eilers
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @export
#' @import grDevices

set_panels = function(rows = 1, cols = 1) {
  # Set number of rows and column for multiple panels
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(rows, cols))
  par(mar = c(4, 3, 2, 2))  # Margins around panels
  par(mgp = c(1.6, 0.6, 0))  # Distance of labels and tick marks
  par(tcl = -0.4)  # Tick length
}

