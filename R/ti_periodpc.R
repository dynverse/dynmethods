#' Inferring trajectories with Periodic PrinCurvee
#'
#' Arguments passed to this function will be used as default parameters for the method.
#'
#' @export
#'
#' @include wrapper_create_description.R
description_periodpc <- create_description(
  name = "Periodic PrinCurve",
  short_name = "periodpc",
  package_loaded = c("stats", "princurve"),
  package_required = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "ndim", default = 3L, lower = 2L, upper = 10L),
    makeIntegerParam(id = "maxit", default = 10L, lower = 0L, upper = 100L)
  ),
  run_fun = "run_periodpc",
  plot_fun = "plot_periodpc"
)

run_periodpc <- function(
  expression,
  ndim = 3,
  maxit = 10
) {
  requireNamespace("stats")
  requireNamespace("princurve")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform PCA dimred
  dimred <- dimred(expression, method = "pca", ndim = ndim)

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
  ) %>% add_cyclic_trajectory(
    pseudotimes = pseudotimes
  ) %>% add_dimred(
    dimred = dimred,
    dimred_trajectory_segments = dimred_trajectory_segments
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_periodpc <- function(prediction) {
  dimred_df <- data.frame(prediction$dimred)
  seg_df <- data.frame(prediction$dimred_trajectory_segments)
  g <- ggplot() +
    geom_point(aes(Comp1, Comp2), dimred_df) +
    geom_segment(aes(x = from_Comp1, xend = to_Comp1, y = from_Comp2, yend = to_Comp2), seg_df, colour = "darkgray")
  process_dynplot(g, prediction$id)
}
