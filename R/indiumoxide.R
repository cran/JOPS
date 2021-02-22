#' An X-ray diffractogram.
#'
#' @docType data
#'
#' @usage data(indiumoxide)
#'
#' @format A matrix with two columns:
#' \describe{
#'   \item{\code{angle}}{the angles (degrees) of diffraction}
#'   \item{\code{count}}{corresponding photon counts.}
#' }
#'
#' @details
#' An X-ray diffractogram of Indium-Tin oxide.
#'
#' These data have been taken from the source of package \code{Diffractometry},
#' which is no longer available from CRAN in binary form.
#'
#' @keywords datasets
#'
#' @source
#' P.L. Davies, U. Gather, M. Meise, D. Mergel, T. Mildenberger (2008).
#' Residual based localization and quantification of peaks in x-ray diffractograms,
#' \emph{Annals of Applied Statistics}, Vol. 2, No. 3, 861-886.
#'

#' @examples
#' angle = indiumoxide[,1]
#' photon = indiumoxide[,2]
#' plot(angle, type = 'l', photon, xlab = 'Angle', ylab = 'Photon count')

"indiumoxide"
