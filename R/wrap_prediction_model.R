#' An abstract data wrapper for TI predictions
#'
#' @inheritParams dynutils::wrap_data
#'
#' @export
wrap_prediction_model <- function(
  cell_ids,
  cell_info = NULL,
  ...
) {
  wrap_data(
    id = random_time_string("TIpred"),
    cell_ids = cell_ids,
    cell_info = cell_info,
    ...
  )
}
