#' Compute the row tensor product of two matrices
#'
#' @description Compute the row tensor product of two matrices with identical numbers of rows.
#'
#' @param X a numeric matrix.
#' @param Y a numeric matrix (if missing, \code{Y = x}).
#'
#' @return The row-wise tensor product of the two matrices.
#'
#' @details
#' The input matrices must have the same number of rows, say \code{m}. If their numbers of columns are \code{n1} and \code{n2},
#' the result is a matrix with \code{m} rows and \code{n1 * n2} columns. Each row of the result is the Kronecker
#' product of the corresponding rows of \code{X} and \code{Y}.
#'
#' @author Paul Eilers
#'
#' @references
#' Eilers, P. H. C. and Currie, I. D. and Durban, M. (2006)
#' Fast and compact smoothing on large multidimensional grids
#' \emph{CSDA} 50, 61--76.
#'
#' @export

rowtens = function(X, Y = X){
  # Row-wise tensor products
  onex = matrix(1, nrow = 1, ncol = ncol(X))
  oney = matrix(1, nrow = 1, ncol = ncol(Y))
  kronecker(X, oney) * kronecker(onex, Y)
}
