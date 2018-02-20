#' Wrap a linear TI prediction
#'
#' This function will generate the milestone_network and progressions.
#'
#' @param pseudotimes The pseudotimes of the \code{cell_ids}
#' @inheritParams wrap_prediction_model
#'
#' @export
wrap_prediction_model_linear <- function(
  cell_ids,
  pseudotimes,
  ...
) {
  .Deprecated("add_linear_trajectory_to_wrapper", package = "dynutils")
}
