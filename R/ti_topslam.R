#' Inferring trajectories with topslam
#'
#' Arguments passed to this function will be used as default parameters for the method.
#'
#' @export
#'
#' @include wrapper_create_description.R
description_topslam <- create_description(
  name = "topslam",
  short_name = "topslam",
  package_loaded = c(),
  package_required = c("topslam"),
  par_set = makeParamSet(
    makeIntegerParam(id = "n_components", lower = 2L, default=2L, upper=10L),
    makeIntegerParam(id = "n_neighbors", lower = 2L, default=10L, upper=100L),
    makeIntegerParam(id = "linear_dims", lower = 0L, default=0L, upper=5L),
    makeNumericParam(id = "max_iters", lower = log(10), default = log(1000), upper = log(10000), trafo = function(x) round(exp(x))),
    makeLogicalVectorParam(id = "dimreds", len = 5, default = rep(TRUE, 5))
  ),
  run_fun = "run_topslam",
  plot_fun = "plot_topslam"
)

#' @importFrom dplyr bind_cols
run_topslam <- function(
  expression,
  start_cells,
  n_components = 2,
  n_neighbors = 10,
  linear_dims = 0,
  max_iters = 200,
  dimreds = rep(TRUE, 5)
) {
  dimreds_vec <- c("t-SNE", "PCA", "Spectral", "Isomap", "ICA")[dimreds]

  start_cell_id <-
    if (!is.null(start_cells)) {
      sample(start_cells, 1)
    } else {
      NULL
    }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run topslam
  out <- topslam::topslam(
    expression = expression,
    start_cell_id = start_cell_id,
    n_components = n_components,
    n_neighbors = n_neighbors,
    linear_dims = linear_dims,
    max_iters = max_iters,
    dimreds = dimreds_vec
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # combine data
  wad <- with(out, bind_cols(wad_grid, wad_energy))
  model <- with(out, bind_cols(data_frame(cell_id = rownames(expression)), space, pseudotime))

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotimes = model$time %>% setNames(rownames(expression)),
    model = model,
    wad = wad
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom viridis scale_colour_viridis
plot_topslam <- function(prediction) {
  ranges <- prediction$model %>% select(starts_with("Comp")) %>% apply(2, range)

  g <- ggplot() +
    # geom_raster(aes(x, y, fill = energy), prediction$wad) +
    # geom_contour(aes(x, y, z = energy, weight = energy), prediction$wad, binwidth = 0.05, color = "black", alpha = 0.4) +
    geom_point(aes(Comp1, Comp2, color = time), prediction$model) +
    scale_fill_gradientn(colors=c("white", "gray20")) +
    viridis::scale_colour_viridis(option = "plasma") +
    labs(colour = "Pseudotime", fill = "Energy") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}
