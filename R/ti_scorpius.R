#' Inferring trajectories with SCORPIUS
#'
#' @inherit ti_angle description
#'
#' @inheritParams SCORPIUS::reduce_dimensionality
#' @inheritParams SCORPIUS::infer_trajectory
#' @param sparse Whether or not to use sparse MDS dimensionality reduction,
#'   for datasets with large amounts of cells.
#' @param distance_method A character string indicating which correlation
#'  coefficient (or covariance) is to be computed. One of "pearson", "kendall", or "spearman".
#'
#' @export
ti_scorpius <- create_ti_method(
  name = "SCORPIUS",
  short_name = "scorpius",
  implementation_id = "scorpius",
  package_loaded = c(),
  package_required = c("SCORPIUS"),
  doi = "10.1101/079509",
  trajectory_types = c("linear", "bifurcation", "convergence", "multifurcation"),
  topology_inference = "fixed",
  type = "algorithm",
  license = "GPL-3",
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
    ),
    list(
      given = "Yvan",
      family = "Saeys",
      email = "yvan.saeys@ugent.be",
      github = "saeyslab"
    )
  ),
  preprint_date = "2016-10-07",
  version = "1.0",
  code_url = "https://github.com/rcannood/SCORPIUS",
  parameters = list(
    distance_method = list(
      type = "discrete",
      default = "spearman",
      values = c("spearman", "pearson", "kendall"),
      description = "A character string indicating which correlation\ncoefficient (or covariance) is to be computed. One of \"pearson\", \"kendall\", or \"spearman\"."),
    ndim = list(
      type = "integer",
      default = 3L,
      upper = 20L,
      lower = 2L,
      description = "The number of dimensions in the new space."
    ),
    k = list(
      type = "integer",
      default = 4L,
      upper = 20L,
      lower = 1L,
      description = "The number of clusters to cluster the data into."
    ),
    thresh = list(
      type = "numeric",
      default = 1e-3,
      upper = 1e5,
      lower = 1e-5,
      distribution = "exponential",
      rate = 1,
      description = "\\code{\\link[princurve]{principal_curve}} parameter: convergence threshhold on shortest distances to the curve"
    ),
    maxit = list(
      type = "integer",
      default = 10L,
      upper = 50L,
      lower = 0L,
      description = "\\code{\\link[princurve]{principal_curve}} parameter: maximum number of iterations"
    ),
    stretch = list(
      type = "numeric",
      default = 0,
      upper = 5,
      lower = 0,
      description = "\\code{\\link[princurve]{principal_curve}} parameter: a factor by which the curve can be extrapolated when points are projected"
    ),
    smoother = list(
      type = "discrete",
      default = "smooth_spline",
      values = c("smooth_spline", "lowess", "periodic_lowess"),
      description = "\\code{\\link[princurve]{principal_curve}} parameter: choice of smoother"
    ),
    sparse = list(
      type = "logical",
      default = FALSE,
      description = "Whether or not to use sparse MDS dimensionality reduction,\nfor datasets with large amounts of cells."
    )
  ),
  run_fun = "dynmethods::run_scorpius",
  plot_fun = "dynmethods::plot_scorpius"
)

run_scorpius <- function(
  expression,
  ndim = 3,
  k = 4,
  distance_method = "spearman",
  thresh = .001,
  maxit = 10,
  stretch = 0,
  smoother = "smooth.spline",
  sparse = FALSE
) {
  requireNamespace("SCORPIUS")

  # if k is too low, turn off clustering
  if (k <= 1) {
    k <- NULL
  }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  space <- SCORPIUS::reduce_dimensionality(
    x = expression,
    dist_fun = function(x, y) dynutils::calculate_distance(x, y, method = distance_method),
    landmark_method = ifelse(sparse, "naive", "none"),
    ndim = ndim,
    num_landmarks = ifelse(nrow(expression) > 1000, 500, nrow(expression))
  )

  # infer a trajectory through the data
  traj <- SCORPIUS::infer_trajectory(
    space,
    k = k,
    thresh = thresh,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother
  )

  # convert trajectory to segments
  dimred_trajectory_segments <-
    cbind(
      traj$path[-nrow(traj$path), , drop = FALSE] %>% magrittr::set_colnames(., paste0("from_comp_", seq_along(colnames(.)))),
      traj$path[-1, , drop = FALSE] %>% magrittr::set_colnames(., paste0("to_comp_", seq_along(colnames(.))))
    )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotime = traj$time
  ) %>% add_dimred(
    dimred = space,
    dimred_trajectory_segments = dimred_trajectory_segments
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom magrittr set_colnames
plot_scorpius <- function(prediction) {
  requireNamespace("SCORPIUS")
  requireNamespace("MASS")
  requireNamespace("RColorBrewer")

  space <- prediction$dimred
  ranges <- apply(space, 2, range)
  maxrange <- apply(ranges, 2, diff) %>% max

  limits <- ranges %>% sweep(1, c(-.5, .5) * maxrange, "+")

  space_df <- space %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    mutate(time = prediction$pseudotime[cell_id])

  seg_df <- prediction$dimred_trajectory_segments %>%
    as.data.frame

  g <- ggplot() +
    geom_path(aes(comp_1, comp_2), alpha = 0, data.frame(ranges %>% sweep(1, c(-.2, .2) * maxrange, "+"))) +
    geom_point(aes(comp_1, comp_2, colour = time), space_df) +
    geom_segment(aes(x = from_comp_1, xend = to_comp_1, y = from_comp_2, yend = to_comp_2), seg_df) +
    scale_colour_gradientn(colours = RColorBrewer::brewer.pal(3, "Dark2")) +
    labs(colour = "Pseudotime") +
    theme(legend.position = c(.92, .12))

  process_dynplot(g, prediction$id, expand = FALSE)
}
