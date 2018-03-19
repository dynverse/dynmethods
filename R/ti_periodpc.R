#' Description for periodpc
#' @export
description_periodpc <- function() create_description(
  name = "Periodic PrinCurve",
  short_name = "periodpc",
  package_loaded = c("princurve"),
  package_required = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "ndim", default = 3L, lower = 2L, upper = 10L),
    makeIntegerParam(id = "maxit", default = 10L, lower = 0L, upper = 100L)
  ),
  properties = c(),
  run_fun = run_periodpc,
  plot_fun = plot_periodpc
)

run_periodpc <- function(
  expression,
  ndim = 3,
  maxit = 10
) {
  requireNamespace("princurve")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform PCA dimred
  dimred <- prcomp(x)$x[,seq_len(ndim)]

  # apply principal curve with periodic lowess smoother
  fit <- princurve::principal.curve(dimred, smoother = "periodic.lowess", maxit = maxit)

  # get pseudotimes
  pseudotimes <- fit$lambda %>% magrittr::set_names(rownames(expression))

  # construct segments
  path <- fit$s[fit$tag, , drop = FALSE]
  dimred_trajectory_segments <- cbind(
    path,
    path[c(seq(2, nrow(path)), 1), ,drop = F]
  ) %>%
    magrittr::set_colnames(c(paste0("from_", colnames(path)), paste0("to_", colnames(path))))

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_cyclic_trajectory_to_wrapper(
    pseudotimes = pseudotimes
  ) %>% add_dimred_to_wrapper(
    dimred = dimred,
    dimred_trajectory_segments = dimred_trajectory_segments
  ) %>% add_timings_to_wrapper(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_periodpc <- function(prediction) {
  dimred_df <- data.frame(prediction$dimred)
  seg_df <- data.frame(prediction$dimred_trajectory_segments)
  g <- ggplot() +
    geom_point(aes(PC1, PC2), dimred_df) +
    geom_segment(aes(x = from_PC1, xend = to_PC1, y = from_PC2, yend = to_PC2), seg_df, colour = "darkgray")
  process_dynplot(g, prediction$id)
}
