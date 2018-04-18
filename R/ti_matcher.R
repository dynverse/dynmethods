#' Description for MATCHER
#' @export
description_matcher <- function() create_description(
  name = "MATCHER",
  short_name = "matcher",
  package_loaded = c(),
  package_required = c("MATCHER"),
  par_set = makeParamSet(
    makeIntegerParam("quantiles", 2, 500, default=50),
    makeDiscreteParam("method", values=c("gp", "linear"), default="linear")
  ),
  properties = c(),
  run_fun = run_matcher,
  plot_fun = plot_matcher
)

#' @import reticulate
run_matcher <- function(
  counts,
  quantiles=50,
  method="gp"
) {
  requireNamespace("MATCHER")

  # load matcher
  use_virtualenv(file.path(find.package("MATCHER"), "venv"))
  pymatcher <- import("pymatcher")

  # setup data
  X <- np_array(np_array(counts))

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run matcher
  m <- pymatcher$matcher$MATCHER(list(X))
  m$infer(quantiles=as.integer(quantiles), method=list("gp"))

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # extract results
  pseudotimes <- m$master_time[[1]][, 1] %>% set_names(rownames(counts))
  sample_master_time <- m$sample_master_time(0L)
  colnames(sample_master_time) <- names(pseudotimes)

  # return output
  prediction <- wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_linear_trajectory_to_wrapper(
    pseudotimes = pseudotimes,
    sample_master_time = sample_master_time
  ) %>% add_timings_to_wrapper(
    tl %>% add_timing_checkpoint("method_afterpostproc")
  )
  prediction
}

plot_matcher <- function(prediction) {
  prediction$sample_master_time %>%
    reshape2::melt(varnames=c("sample_id", "cell_id"), value.name="pseudotime") %>%
    mutate_if(is.factor, as.character) %>%
    left_join(
      tibble(global_pseudotime = prediction$pseudotimes, cell_id = prediction$cell_id),
      "cell_id"
    ) %>%
    ggplot(aes(global_pseudotime, pseudotime)) +
    geom_point(shape="x", size=2) +
    geom_smooth()
}
