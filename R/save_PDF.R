#' Save a plot as a PDF file.
#'
#' @description Save a plot as a PDF file in a (default) folder.
#' The present default is determined by the folder structure for the production of the book.
#'
#' @param fname the file name without the extension PDF (default: \code{scratch}).
#' @param folder the folder for saving PDF plots (default \code{../../Graphs)}.
#' @param show a logical parameter; if TRUE the full file name will be displayed.
#' @param width figure width in inches (default = 6).
#' @param height figure height in inches (default = 4.5).
#' @return save a plot as a PDF file.
#'
#' @author Paul Eilers
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @export
#' @import grDevices
#'
save_PDF = function(fname = "scratch", folder = "../../Graphs", show = T, width = 6, height = 4.5) {
  fpdf = paste(fname, ".pdf", sep = "")
  ffull = file.path(folder, fpdf)
  if (show) cat(ffull, "\n")
  dev.copy2pdf(file = ffull, width = width, height = height)
}

