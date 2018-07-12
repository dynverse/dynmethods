#' Inferring trajectories with Periodic PrinCurvee
#'
#' @inherit ti_angle description
#'
#' @inheritParams dyndimred::dimred
#' @inheritParams princurve::principal.curve
#'
#' @seealso [princurve::principal.curve()]
#'
#' @export
ti_periodpc <- create_ti_method(
  name = "Periodic PrinCurve",
  short_name = "periodpc",
  package_loaded = c("stats", "princurve"),
  package_required = c(),
  trajectory_types = c("cycle", "linear", "bifurcation", "convergence", "multifurcation", "binary_tree", "tree", "acyclic_graph", "graph", "disconnected_graph"),
  topology_inference = "fixed",
  type = "algorithm_test",
  authors = list(
    list(
      given = "Robrecht",
      family = "Cannoodt",
      email = "rcannood@gmail.com",
      ORCID = "0000-0003-3641-729X",
      github = "rcannood"
    ),
    list(
      given = "Wouter",
      family = "Saelens",
      email = "wouter.saelens@ugent.be",
      ORCID = "0000-0002-7114-6248",
      github = "zouter"
    )
  ),
  parameters = list(
    ndim = list(
      type = "integer",
      default = 3,
      lower = 2,
      upper = 10
    ),
    maxit = list(
      type = "integer",
      default = 10,
      lower = 0,
      upper = 100
    )
  ),
  run_fun = "dynmethods::run_periodpc",
  plot_fun = "dynmethods::plot_periodpc"
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
  dimred <- dyndimred::dimred(expression, method = "pca", ndim = ndim)

  # apply principal curve with periodic lowess smoother
  fit <- princurve::principal.curve(dimred, smoother = "periodic.lowess", maxit = maxit)

  # get pseudotime
  pseudotime <- fit$lambda %>% magrittr::set_names(rownames(expression))

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
    pseudotime = pseudotime
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
    geom_point(aes(comp_1, comp_2), dimred_df) +
    geom_segment(aes(x = from_comp_1, xend = to_comp_1, y = from_comp_2, yend = to_comp_2), seg_df, colour = "darkgray")
  process_dynplot(g, prediction$id)
}
