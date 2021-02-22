#' Chromosome G519C18 data
#'
#' @description An extract of the data set \code{G519} in the Bioconductor package \code{Vega}, for chromosome 18.
#'
#' @docType data
#'
#' @usage data(G519C18)
#'
#' @format A dataframe with two columns:
#' \describe{
#'   \item{\code{y}}{Probe position}
#'   \item{\code{x}}{Log R Ratio}. }
#'
#' @keywords datasets
#'
#' @references https://www.bioconductor.org/packages/release/bioc/html/Vega.html
#'
#' @examples
#' plot(G519C18$x, G519C18$y, type = 'l', ylab = 'LRR', xlab = 'Position', main = 'Chromosome 18')

"G519C18"
