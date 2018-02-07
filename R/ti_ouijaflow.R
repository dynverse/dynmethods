#' Description for Ouijaflow
#' @export
description_ouijaflow <- function() create_description(
  name = "ouijaflow",
  short_name = "ouijaflw",
  package_required = c("ouijaflow"),
  package_loaded = c(),
  par_set = makeParamSet(
    makeNumericParam(id = "iter", lower = log(2), default = log(1000), upper = log(50000), trafo = function(x) round(exp(x)))
  ),
  properties = c(),
  run_fun = run_ouijaflow,
  plot_fun = plot_ouija
)

run_ouijaflow <- function(
  expression,
  iter = 1000
) {

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  pseudotimes <- ouijaflow::ouijaflow(expression, iter=iter)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # TIMING: after postproc
  tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

  wrap_prediction_model_linear(
    cell_ids = rownames(expression),
    pseudotimes = pseudotimes
  ) %>% attach_timings_attribute(tl)
}
