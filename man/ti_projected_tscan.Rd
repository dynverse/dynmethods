% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ti_projected_tscan.R
\name{ti_projected_tscan}
\alias{ti_projected_tscan}
\title{Projected TSCAN}
\usage{
ti_projected_tscan(
  minexpr_percent = 0L,
  minexpr_value = 0L,
  cvcutoff = 0L,
  clusternum = c(2L, 9L),
  modelNames = "VVV"
)
}
\arguments{
\item{minexpr_percent}{. Domain: U(0, 1). Default: 0. Format: numeric.}

\item{minexpr_value}{. Domain: U(0, 10). Default: 0. Format: numeric.}

\item{cvcutoff}{. Domain: U(0, 5). Default: 0. Format: numeric.}

\item{clusternum}{. Domain: ( U(2, 20), U(2, 20) ). Default: (2, 9). Format:
integer_range.}

\item{modelNames}{. Domain: {EII, VII, EEI, VEI, EVI, VVI, EEE, EVE, VEE, VVE,
EEV, VEV, EVV, VVV}. Default: VVV. Format: character.}
}
\value{
A TI method wrapper to be used together with
\code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
}
\description{
Will generate a trajectory using \href{https://github.com/zji90/TSCAN}{Projected TSCAN}.

This method was wrapped inside a
\href{https://github.com/dynverse/ti_projected_tscan}{container}.
The original code of this method is available
\href{https://github.com/zji90/TSCAN}{here}.
}
\keyword{method}
