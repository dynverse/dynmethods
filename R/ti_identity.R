#' Description for identity
#'
#' @importFrom dynplot plot_default
#'
#' @export
description_identity <- function() create_description(
  name = "Control: identity",
  short_name = "identity",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeNumericParam(id = "dummy_param", lower = 0, default = 0.5, upper = 1)
  ),
  properties = c(),
  run_fun = run_identity,
  plot_fun = dynplot::plot_default
)

run_identity <- function(counts, task, dummy_param = .5) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = task$cell_ids,
    cell_info = task$cell_info
  ) %>%
    add_trajectory_to_wrapper(
      milestone_ids = task$milestone_ids,
      milestone_network = task$milestone_network,
      divergence_regions = task$divergence_regions,
      progressions = task$progressions
    ) %>%
    add_timings_to_wrapper(
      timings = tl %>% add_timing_checkpoint("method_afterpostproc")
    )
}

