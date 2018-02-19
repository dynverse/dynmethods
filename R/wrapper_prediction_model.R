############### TODO: REMOVE THIS FILE


#' Wrapping TI predictions
#'
#' @param cell_ids the ids of the cells in the trajectory and counts
#' @param milestone_ids the ids of the milestones in the trajectory
#' @param milestone_network the network of the milestones
#' @param milestone_percentages what percentage of milestone is each cell
#' @param progressions what progression does a cell have
#' @param ... extra information to be stored in the wrapper
#'
#' @export
wrap_prediction_model <- function(
  cell_ids,
  milestone_ids,
  milestone_network,
  milestone_percentages = NULL,
  progressions = NULL,
  ...
) {
  .Deprecated("add_trajectory_to_wrapper", package = "dynutils")
}

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
