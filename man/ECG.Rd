% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECG.R
\docType{data}
\name{ECG}
\alias{ECG}
\title{A section of an ECG (electrocardiogram)}
\format{
A data frame with three columns:
\describe{
\item{time}{time in seconds}
\item{resp}{respiration, arbitrary units}
\item{ecg}{ECG, arbitrary units.}
}
}
\source{
https://physionet.org/content/fantasia/1.0.0/
}
\usage{
data(ECG)
}
\description{
The data set includes two signals, respiration and the ECG.
Both signals are distorted by strong 60Hz interference from the mains power.
}
\references{
Iyengar N, Peng C-K, Morin R, Goldberger AL, Lipsitz LA.
Age-related alterations in the fractal scaling of cardiac interbeat interval dynamics.
\emph{Am J Physiol}, 1996; 271: 1078-1084.

Standard citation for PhysioNet:
Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG, Mietus JE, Moody GB, Peng C-K, Stanley HE.
PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals (2003).
Circulation. 101(23):e215-e220.
}
\keyword{datasets}
