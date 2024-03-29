% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ti_gpfates.R
\name{ti_gpfates}
\alias{ti_gpfates}
\title{GPfates}
\usage{
ti_gpfates(
  log_expression_cutoff = 0.5,
  min_cells_expression_cutoff = 0L,
  ndim = 2L
)
}
\arguments{
\item{log_expression_cutoff}{The log expression cutoff. Domain: U(0.5, 5).
Default: 0.5. Format: numeric.}

\item{min_cells_expression_cutoff}{The min expression cutoff. Domain: U(0, 20).
Default: 0. Format: numeric.}

\item{ndim}{Number of dimensions for dimensionality reduction. Domain: U(1, 5).
Default: 2. Format: integer.}
}
\value{
A TI method wrapper to be used together with
\code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
}
\description{
Will generate a trajectory using
\href{https://doi.org/10.1126/sciimmunol.aal2192}{GPfates}.

This method was wrapped inside a
\href{https://github.com/dynverse/ti_gpfates}{container}.
The original code of this method is available
\href{https://github.com/Teichlab/GPfates}{here}.
}
\references{
Lönnberg, T., Svensson, V., James, K.R., Fernandez-Ruiz, D.,
Sebina, I., Montandon, R., Soon, M.S.F., Fogg, L.G., Nair, A.S., Liligeto,
U.N., Stubbington, M.J.T., Ly, L.-H., Bagger, F.O., Zwiessele, M., Lawrence,
N.D., Souza-Fonseca-Guimaraes, F., Bunn, P.T., Engwerda, C.R., Heath, W.R.,
Billker, O., Stegle, O., Haque, A., Teichmann, S.A., 2017. Single-cell RNA-seq
and computational analysis using temporal mixture modeling resolves TH1/TFHfate
bifurcation in malaria. Science Immunology 2, eaal2192.
}
\keyword{method}
