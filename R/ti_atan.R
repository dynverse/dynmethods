#' Description for atan
#' @export
description_atan <- function() create_description(
  name = "Arc-tangent",
  short_name = "atan",
  package_loaded = c(""),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "dimred", default = "pca", values = names(list_dimred_methods()))
  ),
  properties = c(),
  run_fun = run_atan,
  plot_fun = plot_atan
)

run_atan <- function(
  expression,
  dimred = "pca"
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform PCA dimred
  space <- dimred(expression, method = dimred, ndim = 2)

  # transform to pseudotimes using atan2
  pseudotimes <- atan2(dimred[,2], dimred[,1]) / 2 / pi

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
plot_atan <- function(prediction) {
  dimred_df <- data.frame(prediction$dimred, pseudotime = cos(prediction$pseudotimes * 2 * pi))
  g <- ggplot() +
    geom_point(aes(PC1, PC2, colour = pseudotime), dimred_df) +
    viridis::scale_color_viridis()
  process_dynplot(g, prediction$id)
}
