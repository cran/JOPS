% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Complaints.R
\docType{data}
\name{Complaints}
\alias{Complaints}
\title{Environmental complaints from the Rijnomond area of The Netherlands}
\format{
A dataframe with two columns:
\describe{
  \item{\code{freq}}{The daily number of complaints.}
  \item{\code{count}}{The number of days the specific complaint frequency occurred.}
  }
}
\source{
Personal information from Paul Eilers.
}
\usage{
data(Complaints)
}
\description{
Environmental complaints about odors from the Rijnmond region (near Rotterdam in the Netherlands) in 1988.
}
\details{
In 1988, the Rijnmond Environmental Agency registered approximately 20,000
complaints about odors from regional inhabitants.
}
\examples{
plot(Complaints$freq, Complaints$count, type = 'h',
xlab = 'Number of complaints per day', ylab = 'Frequency')
}
\keyword{datasets}
