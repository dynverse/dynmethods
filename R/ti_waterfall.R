#' Inferring trajectories with Waterfall
#'
#' @inherit ti_angle description
#'
#' @param num_clusters Number of clusters
#'
#' @export
ti_waterfll <-  create_ti_method(
  name = "Waterfall",
  short_name = "waterfll", # max 8 chars
  package_loaded = c(),
  package_required = c("Waterfall"),
  par_set = makeParamSet(
    makeIntegerParam(id = "num_clusters", lower = 2L, default = 10L, upper = 20L)
  ),
  run_fun = "dynmethods::run_waterfall",
  plot_fun = "dynmethods::plot_waterfall"
)

run_waterfall <- function(expression, num_clusters = 10) {
  requireNamespace("Waterfall")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run waterfall
  ps <- Waterfall::pseudotimeprog.foo(t(expression), k = num_clusters)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotime = ps$pseudotime %>% setNames(rownames(expression)),
    ps = ps
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_waterfall <- function(prediction) {
  # adapted from Waterfall::plot_waterfall(prediction$ps)
  flow <- attr(prediction$ps, "flow")
  g <- ggplot() +
    geom_point(aes(pseudotime, pseudotime.y, colour = pseudotime), prediction$ps, size = 5) +
    geom_point(aes(PC1, PC2), flow, size = 2, colour = "red") +
    geom_path(aes(PC1, PC2), flow, size = 2, colour = "red") +
    scale_colour_distiller(palette = "RdBu") +
    labs(colour = "Pseudotime") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}
