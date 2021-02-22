#' Inverse link function, used for GLM fitting.
#'
#' @param x scalar, vector, or matrix input.
#' @param link the link function, one of \code{"identity"}, \code{"log"}, \code{"sqrt"},
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"loglog"}, \code{"reciprocal"};
#' quotes are needed (default \code{"identity"}).
#'
#' @return The inverse link function applied to \code{x}.
#' If \code{link} is not in the above list of allowed  names, \code{NULL} will be returned.
#'
#' @export
#'
inverse_link <- function(x, link) {
  switch(link,
         identity = x,
         logit = 1 / (1 + exp(-x)),
         probit = pnorm(x),
         cloglog = (1 - exp(-exp(x))),
         loglog = exp(-exp(-x)),
         sqrt = x ^ 2,
         log = exp(x),
         reciprocal = 1 /x
        )
}

