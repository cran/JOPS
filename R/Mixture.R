#' Mixture Data
#'
#' @description The mixture data were obtained in an unpublished experiment in
#' 2001 by Zhenyu Wang at University of Amsterdam, under the
#' supervision of Age Smilde. We are grateful for the permission
#' to use the data.
#'
#' @details
#' The following instruments and chemicals were used in the
#' experiment: HP 8453 spectrophotometer (Hewlett-Packard, Palo
#' Alto, CA); 2cm closed quartz cuvette with glass thermostatable
#' jacket; Pt-100 temperature sensor; Neslab microprocessor EX-111 circulator bath;
#' UV-visible Chemstation software (Rev A.02.04) on a
#' Hewlett-Packard Vectra XM2 PC.
#'
#' @docType data
#'
#' @usage data(Mixture)
#'
#' @format A list consisting of the following:
#' \describe{
#' \item{\code{fractions}}{a 34 x 3 matrix of mixure fractions (rows sum to unity):
#' \code{Water} (subboiled demi water (self made)),
#' \code{1,2ethanediol} (99.8\% Sigma-Aldrich Germany),
#' \code{3amino1propanol} (99\% Merk Schuchardt Germany)}
#' \item{\code{xspectra}}{spectra array, 34 (observations) x 401 (wavelenths
#' channels) x 12 (temperatures (C):
#' 30, 35, 37.5, 40, 45, 47.5, 50, 55, 60, 62.5, 65, 70 )}
#' \item{\code{wl}}{wavelengths for the spectra, 700 to 1100 (nm), by 1nm.}
#' }
#'
#' @keywords datasets
#'
#' @references  Eilers, P. H. C., and Marx, B. D. (2003). Multivariate calibration with temperature interaction
#' using two-dimensional penalized signal regression. \emph{Chemometrics and
#' Intellegent Laboratory Systems}, 66, 159–174.
#'
#' @references Marx, B. D., Eilers, P. H. C., and Li, B. (2011). Multidimensional single-index signal
#' regression. \emph{Chemometrics and Intelligent Laboratory Systems}, 109(2), 120–130.
#' [see the Appendix within]
#'
#' @references Zhenyou Wang and Age Smilde, Univeristy of Amsterdam,
#' The Netherlands. Personal communication.
"Mixture"
