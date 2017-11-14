#' Description for identity
#' @export
description_identity <- function() create_description(
  name = "identity",
  short_name = "identity",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeNumericParam(id = "dummy_param", lower = 0, default = 0.5, upper = 1)
  ),
  properties = c(),
  run_fun = run_identity,
  plot_fun = plot_default
)

run_identity <- function(counts, task, dummy_param = .5) {
  # return output
  wrap_ti_prediction(
    ti_type = task$ti_type,
    id = "identity",
    cell_ids = task$cell_ids,
    milestone_ids = task$milestone_ids,
    milestone_network = task$milestone_network,
    progressions = task$progressions
  )
}

