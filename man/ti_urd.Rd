% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ti_urd.R
\name{ti_urd}
\alias{ti_urd}
\title{URD}
\usage{
ti_urd(
  knn = 0L,
  sigma.use = 0,
  distance = "euclidean",
  n_floods = 20L,
  stability.div = 10L,
  mp.factor = 1L,
  perplexity = 30L,
  theta = 0.5,
  max_iter = 1000L,
  num.nn = 30L,
  do.jaccard = TRUE,
  optimal.cells.forward = 20L,
  max.cells.back = 40L,
  n.per.tip = 25000L,
  root.visits = 1L,
  max.steps = 25000L,
  n.subsample = 10L,
  divergence.method = "ks",
  cells.per.pseudotime.bin = 80L,
  bins.per.pseudotime.window = 5L,
  p.thresh = 0.01
)
}
\arguments{
\item{knn}{Number of nearest neighbors to use. \code{0} takes a guess. Domain: U(0,
50). Default: 0. Format: integer.}

\item{sigma.use}{Kernel width to use for the diffusion map. \code{0} uses destiny's
global auto-detection procedure. Domain: U(0, 1). Default: 0. Format: numeric.}

\item{distance}{Distance metric to use for determining transition
probabilities. Domain: {euclidean, cosine, rankcor}. Default: euclidean.
Format: character.}

\item{n_floods}{Number of simulations to perform and average. Domain: U(5, 50).
Default: 20. Format: integer.}

\item{stability.div}{Number of simulation subsamplings to calculate. Domain:
U(2, 50). Default: 10. Format: numeric.}

\item{mp.factor}{Retain PCs than are this factor more than the estimated
maximum singular value expected or random data. This is useful in cases when
there are many PCs that have standard deviations just above that expected by
random, which probably represent noise and should be excluded. Domain: U(0,
10). Default: 1. Format: numeric.}

\item{perplexity}{Perplexity parameter for the tSNE. Domain: U(0, 100).
Default: 30. Format: numeric.}

\item{theta}{Speed/accuracy trade-off for Barnes-Hut approximation of tSNE. 0
is exact tSNE, higher is less accurate. Domain: U(0, 1). Default: 0.5. Format:
numeric.}

\item{max_iter}{Number of nearest neighbors to use. \code{0} takes a guess. Domain:
e^U(4.61, 9.21). Default: 1000. Format: integer.}

\item{num.nn}{How many nearest-neighbors to use in the k-nn graph. Domain:
e^U(2.30, 4.61). Default: 30. Format: integer.}

\item{do.jaccard}{Weight edges in the k-nn graph according to their Jaccard
overlap?. Default: TRUE. Format: logical.}

\item{optimal.cells.forward}{The number of cells in the direction specified by
pseudotime.direction at which the logistic should reach 1-asymptote. Domain:
e^U(1.61, 4.61). Default: 20. Format: numeric.}

\item{max.cells.back}{The number of cells in the direction opposite from that
specified by pseudotime.direction at which the logistic should reach asymptote.
Domain: e^U(1.61, 5.30). Default: 40. Format: numeric.}

\item{n.per.tip}{Number of walks to do per tip. Domain: e^U(4.61, 13.82).
Default: 25000. Format: integer.}

\item{root.visits}{Number of steps to take that visit a root.cell before
stopping. Domain: U(1, 5). Default: 1. Format: integer.}

\item{max.steps}{Number of walks to do per tip. Domain: e^U(4.61, 13.82).
Default: 25000. Format: integer.}

\item{n.subsample}{Number of subsamplings to perform for calculating stability.
Domain: e^U(0.69, 4.61). Default: 10. Format: integer.}

\item{divergence.method}{Distance metric to use for determining transition
probabilities. Domain: {ks, preference}. Default: ks. Format: character.}

\item{cells.per.pseudotime.bin}{Approximate number of cells to assign to each
pseudotime bin for branchpoint finding. Domain: e^U(2.30, 6.91). Default: 80.
Format: integer.}

\item{bins.per.pseudotime.window}{Width of moving window in pseudotime used for
branchpoint finding, in terms of bins. Domain: U(2, 20). Default: 5. Format:
integer.}

\item{p.thresh}{P-value threshold to use in determining whether visitation is
significantly different from pairs of tips. Domain: e^U(-11.51, 0.00). Default:
0.01. Format: numeric.}
}
\value{
A TI method wrapper to be used together with
\code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
}
\description{
Will generate a trajectory using
\href{https://doi.org/10.1126/science.aar3131}{URD}.

This method was wrapped inside a
\href{https://github.com/dynverse/ti_urd}{container}.
The original code of this method is available
\href{https://github.com/farrellja/URD}{here}.
}
\references{
Farrell, J.A., Wang, Y., Riesenfeld, S.J., Shekhar, K., Regev, A.,
Schier, A.F., 2018. Single-cell reconstruction of developmental trajectories
during zebrafish embryogenesis. Science 360, eaar3131.
}
\keyword{method}
