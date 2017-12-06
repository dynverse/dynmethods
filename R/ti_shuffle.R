#' Description for shuffled
#'
#' @importFrom dynplot plot_default
#'
#' @export
description_shuffle <- function() create_description(
  name = "shuffle",
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

run_shuffle <- function(counts, task, dummy_param = .5) {
  # permute cell labels
  allcells <- rownames(counts)
  mapper <- setNames(sample(allcells), allcells)
  progressions <- task$progressions %>% mutate(
    cell_id = mapper[cell_id]
  )

  # return output
  wrap_ti_prediction(
    trajectory_type = task$trajectory_type,
    id = "shuffled",
    cell_ids = task$cell_ids,
    milestone_ids = task$milestone_ids,
    milestone_network = task$milestone_network,
    progressions = progressions
  )
}

