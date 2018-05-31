#' Inferring trajectories with MATCHER
#'
#' @inherit ti_angle description
#'
#' @param quantiles How many quantiles to use when computing warp functions (integer)
#' @param method Gaussian process regression or linear interpolation? ("gp" or "linear)
#'
#' @export
ti_matcher <- create_ti_method(
  name = "MATCHER",
  short_name = "matcher",
  package_loaded = c(),
  package_required = c("MATCHER"),
  par_set = makeParamSet(
    makeIntegerParam("quantiles", 2, 500, default = 50),
    makeDiscreteParam("method", values = c("gp", "linear"), default = "linear")
  ),
  run_fun = "dynmethods::run_matcher",
  plot_fun = "dynmethods::plot_matcher"
)

#' @import reticulate
run_matcher <- function(
  counts,
  quantiles = 50,
  method = "gp",
  num_cores = 1
) {
  requireNamespace("MATCHER")
  set_cores(num_cores)

  # load matcher
  use_virtualenv(file.path(find.package("MATCHER"), "venv"))
  pymatcher <- import("pymatcher")

  # setup data
  X <- np_array(np_array(counts))

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run matcher
  m <- pymatcher$matcher$MATCHER(list(X))
  m$infer(quantiles = as.integer(quantiles), method = list(method))

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # extract results
  pseudotime <- m$master_time[[1]][, 1] %>% set_names(rownames(counts))

  # returns "ValueError: A value in x_new is above the interpolation range."
  # sample_master_time <- m$sample_master_time(0L, 1L)
  # colnames(sample_master_time) <- names(pseudotime)

  # return output
  prediction <- wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_linear_trajectory(
    pseudotime = pseudotime#,
    # sample_master_time = sample_master_time
  ) %>% add_timings(
    tl %>% add_timing_checkpoint("method_afterpostproc")
  )
  prediction
}

plot_matcher <- function(prediction) {
  warning("not implemened because sample_master_time is broken")
  plot_default(prediction)
  # prediction$sample_master_time %>%
  #   reshape2::melt(varnames = c("sample_id", "cell_id"), value.name = "pseudotime") %>%
  #   mutate_if(is.factor, as.character) %>%
  #   left_join(
  #     tibble(global_pseudotime = prediction$pseudotime, cell_id = prediction$cell_id),
  #     "cell_id"
  #   ) %>%
  #   ggplot(aes(global_pseudotime, pseudotime)) +
  #   geom_point(shape = "x", size=2) +
  #   geom_smooth()
}
