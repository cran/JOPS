#' Open a graphics window.
#'
#' @description Open a a window for graphics, with specified width and height.
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @param width figure width in inches (default = 6).
#' @param height figure height in inches (default = 4.5).
#' @param kill if TRUE (default) closes all graphics windows. Works only for Windows.
#' @param noRStudioGD if TRUE: do not use the RStudio device (which does not accept width and height).
#' @return open a graphics window.
#' @note Currently only works for Windows!
#'
#' @export
#' @import grDevices

set_window = function(width = 6, height = 4.5, kill = TRUE, noRStudioGD = TRUE) {
  if (kill) graphics.off()
  dev.new(width = width, height = height, noRStudioGD)
  set_panels(1, 1)
}


# To be done. Code found on http://doingbayesiandataanalysis.blogspot.com/2013/01/uniform-r-code-for-opening-saving.html

#------------------------------------------------------------------------
# Modified 2013-Jan-22 8:55pm EST

#openGraph = function( width=7 , height=7 , ... ) {
#  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
#    X11( width=width , height=height , type="cairo" , ... )
#  } else { # Windows OS
#    windows( width=width , height=height , ... )
#  }
#}
