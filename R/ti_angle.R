#' Inferring trajectories with Angle
#'
#' @inherit dynwrap::ti_comp1 description
#'
#' @param dimred A character vector specifying which dimensionality reduction method to use.
#'   See [dyndimred::dimred()] for the list of available dimensionality reduction methods.
#'
#' @export
ti_angle <-
  create_ti_method(
    name = "Angle",
    short_name = "angle",
    package_loaded = c(),
    package_required = c(),
    par_set = makeParamSet(
      makeDiscreteParam(id = "dimred", default = "pca", values = names(dyndimred::list_dimred_methods()))
    ),
    run_fun = "dynmethods::run_angle",
    plot_fun = "dynmethods::plot_angle"
  )

run_angle <- function(
  expression,
  dimred
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform PCA dimred
  space <- dimred(expression, method = dimred, ndim = 2)

  # transform to pseudotime using atan2
  pseudotime <- atan2(space[,2], space[,1]) / 2 / pi + .5

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_cyclic_trajectory(
    pseudotime = pseudotime,
    do_scale_minmax = FALSE
  ) %>% add_dimred(
    dimred = space
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom viridis scale_colour_viridis
plot_angle <- function(prediction) {
  dimred_df <-
    data.frame(prediction$dimred, pseudotime = prediction$pseudotime * 2 * pi)
  g <- ggplot() +
    geom_point(aes(comp_1, comp_2, colour = pseudotime), dimred_df) +
    viridis::scale_color_viridis()
  process_dynplot(g, prediction$id)
}
