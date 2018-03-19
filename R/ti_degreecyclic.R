#' Description for pcaacos
#' @export
description_pcaacos <- function() create_description(
  name = "PCA + Arc cosine",
  short_name = "pcaacos",
  package_loaded = c(""),
  package_required = c(),
  par_set = makeParamSet(
  ),
  properties = c(),
  run_fun = run_pcaacos,
  plot_fun = plot_pcaacos
)

run_pcaacos <- function(
  expression,
  ndim = 3,
  maxit = 10
) {
  requireNamespace("princurve")
  requireNamespace("stats")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform PCA dimred
  dimred <- stats::prcomp(expression)$x[,1:2]

  pseudotimes <- acos(dimred[,1] / sqrt(rowSums(dimred^2))) / (2 * pi)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_cyclic_trajectory_to_wrapper(
    pseudotimes = pseudotimes,
    do_scale_minmax = FALSE
  ) %>% add_dimred_to_wrapper(
    dimred = dimred
  ) %>% add_timings_to_wrapper(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom viridis scale_colour_viridis
plot_pcaacos <- function(prediction) {
  dimred_df <- data.frame(prediction$dimred, pseudotime = cos(prediction$pseudotimes * 2 * pi))
  g <- ggplot() +
    geom_point(aes(PC1, PC2, colour = pseudotime), dimred_df) +
    viridis::scale_color_viridis()
  process_dynplot(g, prediction$id)
}
