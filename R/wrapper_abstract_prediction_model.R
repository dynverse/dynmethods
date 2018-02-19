#' An abstract data wrapper for TI predictions
#'
#' @inheritParams dynutils::abstract_data_wrapper
#'
#' @export
abstract_prediction_model <- function(
  cell_ids,
  cell_info = NULL,
  ...
) {
  abstract_data_wrapper(
    id = random_time_string("TIpred"),
    cell_ids = cell_ids,
    cell_info = cell_info,
    ...
  )
}
