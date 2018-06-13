#' Inferring trajectories with Embeddr
#'
#' @inherit ti_angle description
#'
#' @inheritParams embeddr::embeddr
#' @inheritParams embeddr::fit_pseudotime
#' @inheritParams princurve::principal.curve
#' @param ndim Dimension of the embedded space, default is 2
#' @param nn_pct The percentage of cells to use as tge number of nearest neighbours if kernel == 'nn'.
#'
#' @export
ti_embeddr <- create_ti_method(
  name = "Embeddr",
  short_name = "embeddr",
  package_loaded = c(),
  package_required = c("scaterlegacy", "embeddr"),
  parameters = list(
    ndim = list(
      type = "integer",
      default = 2L,
      upper = 10L,
      lower = 2L,
      description = "Dimension of the embedded space, default is 2"),
    kernel = list(
      type = "discrete",
      default = "nn",
      values = c("nn", "dist", "heat"),
      description = "The choice of kernel. 'nn' will give nearest neighbours, 'dist' gives minimum distance and\n'heat' gives a heat kernel. Discussed in detail in 'Laplacian Eigenmaps and Spectral Techniques for Embedding and Clustering',\nBelkin & Niyogi"),
    metric = list(
      type = "discrete",

      default = "correlation",
      values = c("correlation", "euclidean", "cosine"),
      description = "The metric with which to assess 'closeness' for nearest neighbour selection, one of\n'correlation' (pearson) or 'euclidean'. Default is 'correlation'."),
    nn_pct = list(
      type = "numeric",
      default = 0,
      upper = 1,
      lower = -2,
      description = "The percentage of cells to use as tge number of nearest neighbours if kernel == 'nn'."),
    eps = list(
      type = "numeric",
      default = 0,
      upper = 5,
      lower = -5,
      description = "Maximum distance parameter if kernel == 'dist'"),

    t = list(
      type = "numeric",
      default = 0,
      upper = 5,
      lower = -5,
      description = "'time' for heat kernel if kernel == 'heat'"),
    symmetrize = list(
      type = "discrete",
      default = "mean",
      values = c("mean", "ceil", "floor"),
      description = "How to make the adjacency matrix symmetric. Note that slightly\ncounterintuitively, node i having node j as a nearest neighbour doesn't guarantee node\nj has node i. There are several ways to get round this:\n\\itemize{\n\\item \\code{mean} If the above case occurs make the link weight 0.5 so the adjacency matrix becomes \\eqn{0.5(A + A')}\n\\item \\code{ceil} If the above case occurs set the link weight to 1 (ie take the ceiling of the mean case)\n\\item \\code{floor} If the above case occurs set the link weight to 0 (ie take the floor of the mean case)\n}"),

    measure_type = list(
      type = "discrete",
      default = "unorm",
      values = c("unorm", "norm"),
      description = "Type of laplacian eigenmap, which corresponds to the constraint on the eigenvalue problem. If\ntype is 'unorm' (default), then the graph measure used is the identity matrix, while if type is 'norm' then the measure\nused is the degree matrix."),
    thresh = list(
      type = "numeric",
      default = 1e-3,
      upper = 1e5,
      lower = 1e-5,
      distribution = "exponential",
      rate = 1,
      description = "convergence threshold on shortest distances to the curve."),
    maxit = list(

      type = "integer",
      default = 10L,
      upper = 50L,
      lower = 0L,
      description = "maximum number of iterations."),
    stretch = list(
      type = "numeric",
      default = 2,
      upper = 5,
      lower = 0,
      description = "a factor by which the curve can be extrapolated when\npoints are projected.  Default is 2 (times the last segment\nlength). The default is 0 for \\code{smoother} equal to\n\\code{\"periodic_lowess\"}."),
    smoother = list(
      type = "discrete",
      default = "smooth.spline",
      values = c("smooth.spline", "lowess",
                 "periodic.lowess"),
      description = "choice of smoother. The default is\n\\code{\"smooth_spline\"}, and other choices are \\code{\"lowess\"} and\n\\code{\"periodic_lowess\"}. The latter allows one to fit closed curves.\nBeware, you may want to use \\code{
iter = 0} with \\code{lowess()}.")
  ),
  run_fun = "dynmethods::run_embeddr",
  plot_fun = "dynmethods::plot_embeddr"
)

run_embeddr <- function(
  counts,
  ndim,
  kernel,
  metric,
  nn_pct,
  eps,
  t,
  symmetrize,
  measure_type,
  thresh,
  maxit,
  stretch,
  smoother
) {
  requireNamespace("scaterlegacy")
  requireNamespace("embeddr")

  # calculate nn param
  nn <- max(round(log(nrow(counts)) * nn_pct), 9)

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # load data in scaterlegacy
  sce <- scaterlegacy::newSCESet(countData = t(counts))

  # run embeddr
  sce <- embeddr::embeddr(
    sce,
    kernel = kernel,
    metric = metric,
    nn = nn,
    eps = eps,
    t = t,
    symmetrize = symmetrize,
    measure_type = measure_type,
    p = ndim
  )

  # fit pseudotime
  sce <- embeddr::fit_pseudotime(
    sce,
    thresh = thresh,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # construct milestone network
  pseudotime <- as(sce@phenoData, "data.frame")$pseudotime %>%
    setNames(rownames(counts))

  # creating extra output for visualisation purposes
  dimred_cells <- sce@reducedDimension

  traj <- as(sce@phenoData, "data.frame") %>%
    dplyr::arrange(pseudotime) %>%
    select(starts_with("trajectory_")) %>%
    as.matrix()

  dimred_trajectory_segments <- cbind(
    traj[-nrow(traj), , drop = F],
    traj[-1, , drop = F]
  )
  colnames(dimred_trajectory_segments) <- c(
    paste0("from_", colnames(dimred_cells)),
    paste0("to_", colnames(dimred_cells))
  )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_linear_trajectory(
    pseudotime = pseudotime
  ) %>% add_dimred(
    dimred = dimred_cells,
    dimred_trajectory_segments = dimred_trajectory_segments
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_embeddr <- function(prediction) {
  sample_df <- prediction$dimred %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    mutate(time = prediction$pseudotime)
  traj_df <- prediction$dimred_trajectory_segments %>%
    as.data.frame()
  g <- ggplot() +
    geom_point(aes(component_1, component_2, fill = time), sample_df, pch = 21, alpha = .65, size = 3.5) +
    geom_segment(aes(x = from_component_1, xend = to_component_1, y = from_component_2, yend = to_component_2),
                 traj_df, size = 1.5, alpha = 0.8, linetype = 2) +
    scale_fill_distiller(palette = "YlOrRd") +
    scale_colour_distiller(palette = "YlOrRd") +
    theme(legend.position = "none")
  process_dynplot(g, prediction$id)
}
