#' Compute a second order circular differencing matrix
#'
#' @description Compute difference matrix used for circular penalities.
#'
#' @param n number of rows (and columns) of the square differencing matrix.
#'
#' @return A square matrix with \code{n} rows and columns.
#'
#' @author Paul Eilers
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, \emph{SORT}, 39(2): 149-186.
#'
#' @examples
#' # Compare standard and circular differencing matrix
#' n = 8
#' D1 = diff(diag(n), diff = 2)
#' D2 = cdiff(n)
#' oldpar = par(no.readonly = TRUE)
#' on.exit(par(oldpar))
#' par(mfrow = c(1, 2))
#' image(t(D1))
#' title('Linear differencing matrix')
#' image(t(D2))
#' title('Circular differencing matrix')
#'
#'
#' @export

cdiff = function(n) {
  # Compute cyclic difference matrix
  D2 = matrix(0, n, n + 2)
  p = c(1, -2 * cos(2 * pi/n), 1)
  for (k in 1:n) D2[k, (0:2) + k] = p
  D = D2[, 2:(n + 1)]
  D[, 1] = D[, 1] + D2[, n + 2]
  D[, n] = D[, n] + D2[, 1]
  D
}
