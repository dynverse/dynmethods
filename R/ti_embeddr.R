#' Inferring trajectories with Embeddr
#'
#' Arguments passed to this function will be used as default parameters for the method.
#'
#' @export
#'
#' @include wrapper_create_description.R
description_embeddr <- create_description(
  name = "Embeddr",
  short_name = "embeddr",
  package_loaded = c(),
  package_required = c("scaterlegacy", "embeddr"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "kernel", default = "nn", values = c("nn", "dist", "heat")),
    makeDiscreteParam(id = "metric", default = "correlation", values = c("correlation", "euclidean", "cosine")),
    makeNumericParam(id = "nn_pct", lower = -2, upper = log10(10), default = 0, trafo = function(x) 10^x),
    makeNumericParam(id = "eps", lower = -5, upper = 5, default = 0, trafo = function(x) 10^x),
    makeNumericParam(id = "t", lower = -5, upper = 5, default = 0, trafo = function(x) 10^x),
    makeDiscreteParam(id = "symmetrize", default = "mean", values = c("mean", "ceil", "floor")),
    makeDiscreteParam(id = "measure_type", default = "unorm", values = c("unorm", "norm")),
    makeIntegerParam(id = "p", lower = 2L, upper = 10L, default = 2L),
    makeNumericParam(id = "thresh", lower = -5, upper = 5, default = -3, trafo = function(x) 10^x),
    makeIntegerParam(id = "maxit", lower = 0L, upper = 50L, default = 10L),
    makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 2),
    makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "lowess", "periodic.lowess"))
  ),
  run_fun = "run_embeddr",
  plot_fun = "plot_embeddr"
)

run_embeddr <- function(
  counts,
  kernel = "nn",
  metric = "correlation",
  nn_pct = 1,
  eps = 1,
  t = 1,
  symmetrize = "mean",
  measure_type = "unorm",
  p = 2,
  thresh = .001,
  maxit = 10,
  stretch = 2,
  smoother = "smooth.spline"
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
    p = p
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
  pseudotimes <- as(sce@phenoData, "data.frame")$pseudotime %>%
    setNames(rownames(counts))

  # creating extra output for visualisation purposes
  dimred_cells <- sce@reducedDimension

  traj <- as(sce@phenoData, "data.frame") %>%
    arrange(pseudotime) %>%
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
    pseudotimes = pseudotimes
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
    mutate(time = prediction$pseudotimes)
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
