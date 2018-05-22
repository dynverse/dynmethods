#' Description for reCAT
#' @export
description_recat <- function() create_description(
  name = "recat",
  short_name = "recat",
  package_loaded = c(),
  package_required = c("reCAT"),
  par_set = makeParamSet(
    makeIntegerParam(id = "TSPFold", default = 2, lower=2, upper=10),
    makeIntegerParam(id = "beginNum", default = 10, lower=2, upper=20),
    makeIntegerParam(id = "endNum", default = 15, lower=2, upper=20),
    makeIntegerParam(id = "step_size", default = 2, lower=2, upper=20),
    makeIntegerParam(id = "base_cycle_range_start", default = 6, lower=2, upper=20),
    makeIntegerParam(id = "base_cycle_range_end", default = 9, lower=2, upper=20),
    makeIntegerParam(id = "max_num", default = 300, lower=100, upper=500),
    makeDiscreteParam(id = "clustMethod", default = "GMM", values = c("GMM", "Pam", "Kmeans"))
  ),
  properties = c(),
  run_fun = run_recat,
  plot_fun = plot_recat
)

run_recat <- function(
  expression,
  num_cores = 1,
  TSPFold = 2,
  beginNum = 10,
  endNum = 15,
  step_size = 2,
  base_cycle_range_start = 6,
  base_cycle_range_end = 9,
  max_num = 300,
  clustMethod = "GMM"
) {
  requireNamespace("reCAT")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run reCAT
  result <- reCAT::bestEnsembleComplexTSP(
    test_exp = expression,
    TSPFold = TSPFold,
    beginNum = beginNum,
    endNum = endNum,
    base_cycle_range = base_cycle_range_start:base_cycle_range_end,
    step_size = step_size,
    max_num = max_num,
    clustMethod = clustMethod,
    threads = num_cores,
    output = FALSE
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  pseudotimes <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ]

  # wrap
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_cyclic_trajectory(
    pseudotimes = pseudotimes
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_recat <- function(prediction) {
  # TODO
}
