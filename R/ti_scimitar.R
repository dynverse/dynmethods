#' Inferring trajectories with SCIMITAR
#'
#' @inherit ti_identity description
#'
#' @export
#'
#' @include wrapper_create_ti_method.R
ti_scimitar <- create_ti_method(
  name = "SCIMITAR",
  short_name = "scimitar",
  package_loaded = c(),
  package_required = c("SCIMITAR"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "covariance_type", values = c('diag', 'spherical', 'full'), default = "diag"),
    makeIntegerParam(id = "degree", lower = 1, upper = 20, default = 3),
    makeNumericParam(id = "step_size", lower = 0.01, upper = 0.1, default = 0.07),
    makeDiscreteParam(id = "cov_estimator", values = c("identity", "diag", "sample", "global", "glasso", "corpcor", "average"), default = "corpcor"),
    makeNumericParam(id = "cov_reg", lower = 0.01, upper = 0.1, default = 0.05),
    makeIntegerParam(id = "max_iter", lower = 1, upper = 20, default = 3)
  ),
  run_fun = "run_scimitar",
  plot_fun = "plot_scimitar"
)

#' @importFrom readr read_csv
#' @importFrom utils write.table
run_scimitar <- function(
  expression,
  covariance_type = "diag",
  degree = 3,
  step_size = 0.07,
  cov_estimator = "corpcor",
  cov_reg = 0.05,
  max_iter = 3,
  num_cores = 1,
  verbose = FALSE
) {
  requireNamespace("SCIMITAR")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  pseudotime <- SCIMITAR::SCIMITAR(
    expression = expression,
    covariance_type = covariance_type,
    degree = degree,
    step_size = step_size,
    cov_estimator = cov_estimator,
    cov_reg = cov_reg,
    max_iter = max_iter,
    num_cores = num_cores,
    verbose = verbose
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # get percentages of end milestones by multiplying the phi with the pseudotime
  progressions <- pseudotime %>%
    mutate(
      from = "M1",
      to = "M2",
      percentage = dynutils::scale_minmax(pseudotime)
    ) %>%
    select(-pseudotime)

  #  create milestone network
  milestone_network <- data_frame(
    from = "M1",
    to = "M2",
    length = 1,
    directed = TRUE
  )
  milestone_ids <- c("M1", "M2")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    divergence_regions = NULL,
    pseudotime = pseudotime
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_scimitar <- function(prediction) {
  stop("TODO...")
}
