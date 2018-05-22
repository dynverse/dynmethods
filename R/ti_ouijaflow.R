#' Description for Ouijaflow
#' @export
description_ouijaflw <- function() create_description(
  name = "ouijaflow",
  short_name = "ouijaflw",
  package_required = c("ouijaflow"),
  package_loaded = c(),
  par_set = makeParamSet(
    makeNumericParam(id = "iter", lower = log(2), default = log(1000), upper = log(50000), trafo = function(x) round(exp(x)))
  ),
  properties = c(),
  run_fun = run_ouijaflow,
  plot_fun = plot_ouijaflow
)

run_ouijaflow <- function(
  expression,
  iter = 1000
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run ouijaflow
  pseudotimes <- ouijaflow::ouijaflow(expression, iter = iter)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotimes = pseudotimes %>% setNames(rownames(expression))
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_ouijaflow <- function(prediction) {
  # TODO
}
