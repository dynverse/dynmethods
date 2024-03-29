% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ti_celltree_maptpx.R
\name{ti_celltree_maptpx}
\alias{ti_celltree_maptpx}
\title{cellTree maptpx}
\usage{
ti_celltree_maptpx(
  method = "maptpx",
  sd_filter = 0.5,
  absolute_width = 0L,
  width_scale_factor = 1.5,
  outlier_tolerance_factor = 0.1,
  rooting_method = "null",
  num_topics = c(2L, 15L),
  tot_iter = 200L,
  tolerance = 0.05
)
}
\arguments{
\item{method}{LDA inference method to use. Domain: {maptpx}. Default: maptpx.
Format: character.}

\item{sd_filter}{Standard-deviation threshold below which genes should be
removed from the data. Domain: e^U(-4.61, 1.61). Default: 0.5. Format: numeric.}

\item{absolute_width}{Distance threshold below which a cell vertex is
considered to be attached to a backbone vertex (see paper for more details). By
default, this threshold is computed dynamically, based on the distance
distribution for each branch. Domain: U(0, 1). Default: 0. Format: numeric.}

\item{width_scale_factor}{A scaling factor for the dynamically-computed
distance threshold (ignored if absolute_width is provided). Higher values will
result in less branches in the backbone tree, while lower values might lead to
a large number of backbone branches. Domain: e^U(-2.30, 4.61). Default: 1.5.
Format: numeric.}

\item{outlier_tolerance_factor}{Proportion of vertices, out of the total number
of vertices divided by the total number of branches, that can be left at the
end of the backbone tree-building algorithm. Domain: e^U(-9.21, 6.91). Default:
0.1. Format: numeric.}

\item{rooting_method}{Method used to root the backbone tree. Must be either
NULL or one of ‘longest.path’, ‘center.start.group’ or ‘average.start.group’.
‘longest.path’ picks one end of the longest shortest-path between two vertices.
’center.start.group’ picks the vertex in the starting group with lowest
mean-square-distance to the others. ‘average.start.group’ creates a new
artificial vertex, as the average of all cells in the starting group. If no
value is provided, the best method is picked based on the type of grouping and
start group information available. Domain: {longest.path, center.start.group,
average.start.group, null}. Default: null. Format: character.}

\item{num_topics}{Number of topics to fit in the model. Domain: ( U(2, 10),
U(10, 30) ). Default: (2, 15). Format: integer_range.}

\item{tot_iter}{Numeric parameters (optional) forwarded to the chosen LDA
inference method's contol class. Domain: e^U(9.21, 13.82). Default: 200.
Format: numeric.}

\item{tolerance}{Numeric parameters (optional) forwarded to the chosen LDA
inference method's contol class. Domain: e^U(-6.91, -0.69). Default: 0.05.
Format: numeric.}
}
\value{
A TI method wrapper to be used together with
\code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
}
\description{
Will generate a trajectory using \href{https://doi.org/10.1186/s12859-016-1175-6}{cellTree maptpx}.

This method was wrapped inside a
\href{https://github.com/dynverse/ti_celltree_maptpx}{container}.
The original code of this method is available
\href{https://bioconductor.org/packages/release/bioc/html/cellTree.html}{here}.
}
\references{
duVerle, D.A., Yotsukura, S., Nomura, S., Aburatani, H., Tsuda, K.,
2016. CellTree: an R/bioconductor package to infer the hierarchical structure
of cell populations from single-cell RNA-seq data. BMC Bioinformatics 17.
}
\keyword{method}
