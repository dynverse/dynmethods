#' Inferring trajectories with \code{angle}
#'
#' Arguments passed to this function will be used as default parameters for the method.
#'
#' @param dimred A character vector specifying which dimensionality reduction method to use.
#'   See \code{\link{list_dimred_methods}} for the list of available dimensionality reduction methods.
#'
#' @export
description_angle <- function(
  dimred = "pca"
) {
  create_description(
    name = "Angle",
    short_name = "angle",
    package_loaded = c(),
    package_required = c(),
    par_set = makeParamSet(
      makeDiscreteParam(id = "dimred", default = dimred, values = names(list_dimred_methods()))
    ),
    properties = c(),
    run_fun = run_angle,
    plot_fun = plot_angle
  )
}
run_angle <- function(
  expression,
  dimred
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform PCA dimred
  space <- dimred(expression, method = dimred, ndim = 2)

  # transform to pseudotimes using atan2
  pseudotimes <- atan2(space[,2], space[,1]) / 2 / pi + .5

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_cyclic_trajectory(
    pseudotimes = pseudotimes,
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
    data.frame(prediction$dimred, pseudotime = prediction$pseudotimes * 2 * pi)
  g <- ggplot() +
    geom_point(aes(Comp1, Comp2, colour = pseudotime), dimred_df) +
    viridis::scale_color_viridis()
  process_dynplot(g, prediction$id)
}
