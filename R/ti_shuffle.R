#' Description for shuffled
#'
#' @importFrom dynplot plot_default
#'
#' @export
description_shuffle <- function() create_description(
  name = "Control: shuffle",
  short_name = "shuffle",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeNumericParam(id = "dummy_param", lower = 0, default = 0.5, upper = 1)
  ),
  properties = c(),
  run_fun = run_shuffle,
  plot_fun = dynplot::plot_default
)

run_shuffle <- function(
  counts,
  task,
  dummy_param = .5
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # permute cell labels
  allcells <- rownames(counts)
  mapper <- setNames(sample(allcells), allcells)
  progressions <- task$progressions %>% mutate(
    cell_id = mapper[cell_id]
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = task$cell_ids
  ) %>% add_trajectory_to_wrapper(
    milestone_ids = task$milestone_ids,
    milestone_network = task$milestone_network,
    progressions = progressions,
    divergence_regions = task$divergence_regions
  ) %>% add_timings_to_wrapper(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

