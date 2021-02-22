#' Environmental complaints from the Rijnomond area of The Netherlands
#'
#' @description Environmental complaints about odors from the Rijnmond region (near Rotterdam in the Netherlands) in 1988.
#'
#' @docType data
#'
#' @usage data(Complaints)
#'
#' @format A dataframe with two columns:
#' \describe{
#'   \item{\code{freq}}{The daily number of complaints.}
#'   \item{\code{count}}{The number of days the specific complaint frequency occurred.}
#'   }
#'
#' @details In 1988, the Rijnmond Environmental Agency registered approximately 20,000
#' complaints about odors from regional inhabitants.
#'
#' @keywords datasets
#'
#' @source
#' Personal information from Paul Eilers.
#'
#' @examples
#' plot(Complaints$freq, Complaints$count, type = 'h',
#' xlab = 'Number of complaints per day', ylab = 'Frequency')

"Complaints"
