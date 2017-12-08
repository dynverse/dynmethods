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
  trajectory_type,
  cell_ids,
  milestone_ids,
  milestone_network,
  milestone_percentages = NULL,
  progressions = NULL,
  ...
) {
  abstract_data_wrapper(
    type = "ti_pred",
    id = random_time_string("TIpred"),
    cell_ids = cell_ids,
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    progressions = progressions,
    ...
  )
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
  pseudotimes <- scale_minmax(pseudotimes)
  milestone_ids <- c("milestone_A", "milestone_B")
  milestone_network <- data_frame(
    from = milestone_ids[[1]],
    to = milestone_ids[[2]],
    length = 1,
    directed = FALSE
  )
  progressions <- data_frame(
    cell_id = cell_ids,
    from = milestone_ids[[1]],
    to = milestone_ids[[2]],
    percentage = pseudotimes
  )

  wrap_prediction_model(
    cell_ids = cell_ids,
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    pseudotimes = pseudotimes,
    ...
  )
}
