#' Inferring trajectories with Control: identity
#'
#' This control method will return the gold standard.
#'
#' @param dummy_param This parameter does not do anything.
#'
#' @export
#'
#' @include wrapper_create_description.R
description_identity <- create_description(
  name = "Control: identity",
  short_name = "identity",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeNumericParam(id = "dummy_param", lower = 0, default = 0.5, upper = 1)
  ),
  run_fun = "run_identity",
  plot_fun = "dynplot::plot_default"
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
  ) %>% add_trajectory(
    milestone_ids = task$milestone_ids,
    milestone_network = task$milestone_network,
    divergence_regions = task$divergence_regions,
    progressions = task$progressions
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

