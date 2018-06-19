#' Inferring trajectories with Waterfall
#'
#' @inherit ti_angle description
#'
#' @param num_clusters Number of clusters
#'
#' @export
ti_waterfall <-  create_ti_method(
  name = "Waterfall",
  short_name = "waterfall",
  package_loaded = c(),
  package_required = c("Waterfall"),
  doi = "10.1016/j.stem.2015.07.013",
  trajectory_types = "linear",
  topology_inference = "fixed",
  type = "algorithm",
  authors = list(
    list(
      given = "Jaehoon",
      family = "Shin",
      email = "shin@jhmi.edu"
    ),
    list(
      given = "Hongjun",
      family = "Song",
      email = "shongju1@jhmi.edu"
    )
  ),
  publication_date = "2015-09-03",
  code_url = "http://www.cell.com/cms/attachment/2038326541/2052521637/mmc9.zip",
  parameters = list(
    num_clusters = list(type = "integer", lower = 2, default = 10, upper = 20)
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
