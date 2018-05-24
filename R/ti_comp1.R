#' Inferring trajectories with Component 1
#'
#' @inherit ti_angle description
#'
#' @param dimred A character vector specifying which dimensionality reduction method to use.
#'   See \code{\link{list_dimred_methods}} for the list of available dimensionality reduction methods.
#'
#' @export
#'
#' @include wrapper_create_ti_method.R
ti_compone <- create_ti_method(
  name = "Component 1",
  short_name = "comp1",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "dimred", default = "pca", values = names(list_dimred_methods())),
    makeIntegerParam(id="ndim", default = 2, lower=2, upper=30)
  ),
  properties = c(),
  run_fun = "run_compone",
  plot_fun = "plot_compone"
)

run_compone <- function(
  expression,
  ndim,
  dimred
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  space <- dimred(expression, method = dimred, ndim = ndim)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotimes = space[,1] %>% setNames(rownames(expression))
  ) %>% add_dimred(
    dimred = space
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom viridis scale_colour_viridis
plot_compone <- function(prediction) {
  g <- ggplot() +
    geom_point(aes(Comp1, Comp2, color = Comp1), data.frame(prediction$dimred)) +
    viridis::scale_colour_viridis(option = "plasma") +
    labs(colour = "Pseudotime") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}
