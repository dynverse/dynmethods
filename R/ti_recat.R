#' Inferring trajectories with reCAT
#'
#' @inherit ti_angle description
#'
#' @param TSPFold No documentation provided by authors
#' @param beginNum No documentation provided by authors
#' @param endNum No documentation provided by authors
#' @param step_size Determines the number of k to skip in your consensus path, ie if step_size = 2, then reCAT would only calculate and merge the paths for k = 12, 14, 16, 18, …, n-2, n. We recommend step_size of up to a maximum of 5 while preserving the performance of reCAT. Usually a step_size of 2 (by default) would suffice and bigger steps are recommended for larger datasets (>1000 cells) in order to reduce computational time.
#' @param base_cycle_range_start The minimal number of four k’s for computing the reference cycle mentioned in the manuscript. Can be set to 6 or 7
#' @param base_cycle_range_end The maximal number of four k’s for computing the reference cycle mentioned in the manuscript. Can be set to 6 or 7
#' @param max_num No documentation provided by authors
#' @param clustMethod No documentation provided by authors
#'
#' @export
#'
#' @include wrapper_create_ti_method.R
ti_recat <- create_ti_method(
  name = "reCAT",
  short_name = "recat",
  package_loaded = c(),
  package_required = c("reCAT"),
  par_set = makeParamSet(
    makeIntegerParam(id = "TSPFold", default = 2, lower=2, upper=10),
    makeIntegerParam(id = "beginNum", default = 10, lower=2, upper=20),
    makeIntegerParam(id = "endNum", default = 15, lower=2, upper=20),
    makeIntegerParam(id = "step_size", default = 2, lower=2, upper=20),
    makeIntegerParam(id = "base_cycle_range_start", default = 6, lower=6, upper=7),
    makeIntegerParam(id = "base_cycle_range_end", default = 9, lower=9, upper=10),
    makeIntegerParam(id = "max_num", default = 300, lower=100, upper=500),
    makeDiscreteParam(id = "clustMethod", default = "GMM", values = c("GMM", "Pam", "Kmeans"))
  ),
  run_fun = "run_recat",
  plot_fun = "plot_recat"
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

  pseudotimes <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ] %>% set_names(rownames(expression))

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
