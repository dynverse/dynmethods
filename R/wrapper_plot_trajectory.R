#' Visualise a trajectory by the visualisation specific for the method that created it
#'
#' @param prediction The trajectory model to plot
#' @param method The method that generated the trajectory.
#' @param ... Extra parameters for the plotting function.
#' @export
plot_trajectory <- function(
  prediction,
  method,
  ...
) {
  testthat::expect_true(is_prediction(prediction))
  testthat::expect_true(is_ti_method(method))

  method$plot_fun(prediction, ...)
}
