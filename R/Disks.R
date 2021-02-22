#' Prices of hard disk drives
#'
#' @description Prices and capacities of hard disk drives,
#' as advertised in a Dutch computer monthly in 1999.
#' Prices are given in Dutch guilders; the Euro did not yet exist.
#'
#' @docType data
#'
#' @usage data(Disks)
#'
#' @format A dataframe with six columns:
#' \describe{
#'   \item{\code{Year}}{1999-2000}
#'   \item{\code{Month}}{month, 1-12}
#'   \item{\code{Size}}{capacity in Gb}
#'   \item{\code{Buffer}}{buffer size (Mb)}
#'   \item{\code{RPM}}{rotating speed (rpm)}
#'   \item{\code{PriceDG}}{in Dutch Guilders, divide by 2.2 for Euro.}}
#'
#' @keywords datasets
#'
#' @source
#' Personal information from Paul Eilers.
#'

"Disks"
