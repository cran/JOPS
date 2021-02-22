#' Create a matrix of counts.
#'
#' @description Count the number of occurrences of pairs of positive integers in two vectors,
#' producing a matrix.
#'
#' @param xb a vector of integers.
#' @param yb a vector of integers.
#' @param nb a vector of length 2 that provides the number of bins for the 2D histogram on \code{x} and \code{y}.
#'
#' @return A matrix with \code{nb[1]} rows and \code{nb[2]} columns with counts.
#' It serves as the input for two-dimensional histogram smoothing.
#'
#' @details This function builds a two-dimensional histogram, based on two two vectors of bin numbers (obtained
#' with \code{binit}). Rows where \code{x[i] > nb[1]}  or \code{y[i] > nb[2]} are discarded without a warning.
#'
#' @export

count2d <- function(xb, yb, nb) {
  H <- matrix(0, nb[1], nb[2])
  sel = which(xb <= nb[1] & yb <= nb[2])
  for (i in 1:length(sel)) {
    j = sel[i]
    H[xb[j], yb[j]] = H[xb[j], yb[j]] + 1
  }
  H
}
