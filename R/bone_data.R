#' Spinal bone relative mineral density
#'
#' @description Relative spinal bone mineral density measurements on 261 North
#' American adolescents.  Each value is the difference in \code{spnbmd}
#' taken on two consecutive visits, divided by the average. The age is
#' the average age over the two visits.
#'
#' @docType data
#'
#' @usage data(bone_data)
#'
#' @format A dataframe with four columns:
#' \describe{
#'   \item{\code{idnum}}{ID of the child}
#'   \item{\code{age}}{age}
#'   \item{\code{gender}}{male or female}
#'   \item{\code{spnbmd}}{Relative Spinal bone mineral density.}
#' }
#'

#' @keywords datasets
#'
#' @references Bachrach, L.K., Hastie, T., Wang, M.-C., Narasimhan, B., Marcus, R. (1999). Bone Mineral
#' Acquisition in Healthy Asian, Hispanic, Black and Caucasian Youth. A
#' Longitudinal Study. \emph{J Clin Endocrinol Metab} 84, 4702-12.
#'
#' @references Eilers, P.H.C. and Marx, B.D. (2021). \emph{Practical Smoothing, The Joys of
#' P-splines.} Cambridge University Press.
#'
#' @source  https://web.stanford.edu/~hastie/ElemStatLearn/datasets/bone.data
#'

"bone_data"
