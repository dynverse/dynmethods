#' Inferring trajectories with SCORPIUS
#'
#' Arguments passed to this function will be used as default parameters for the method.
#'
#' @inheritParams SCORPIUS::correlation_distance
#' @inheritParams SCORPIUS::reduce_dimensionality
#' @inheritParams SCORPIUS::infer_trajectory
#' @param sparse Whether or not to use sparse MDS dimensionality reduction,
#'   for datasets with large amounts of cells.
#'
#' @rdname scorpius
#'
#' @include wrapper_create_description.R
abstract_scorpius_description <- function(short_name) {
  name <- c(
    "scorpius" = "SCORPIUS",
    "scorspar" = "SCORPIUS sparse"
  )[short_name] %>% setNames(NULL)

  create_description(
    name = name,
    short_name = short_name,
    package_loaded = c(),
    package_required = c("SCORPIUS"),
    par_set = makeParamSet(
      makeDiscreteParam(id = "distance_method", default = "spearman", values = c("spearman", "pearson", "kendall")),
      makeIntegerParam(id = "ndim", lower = 2L, default = 3L, upper = 20L),
      makeIntegerParam(id = "k", lower = 1L, default = 4L, upper = 20L),
      makeNumericParam(id = "thresh", lower = -5, upper = 5, default = -3, trafo = function(x) 10^x),
      makeIntegerParam(id = "maxit", lower = 0L, upper = 50L, default = 10L),
      makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 0),
      makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "lowess", "periodic.lowess")),
      makeLogicalParam(id = "sparse", default = short_name == "scorspar")
    ),
    run_fun = "run_scorpius",
    plot_fun = "plot_scorpius"
  )
}


#' @rdname scorpius
#' @export
description_scorpius <- abstract_scorpius_description("scorpius")

#' @rdname scorpius
#' @export
description_scorpius_sparse <- abstract_scorpius_description("scorspar")

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

  if (!sparse) {
    # calculate distances between cells
    dist <- SCORPIUS::correlation_distance(expression, method = distance_method)

    # perform dimensionality reduction
    space <- SCORPIUS::reduce_dimensionality(dist, ndim = ndim)
  } else {
    dist_fun <- function(x, y) SCORPIUS::correlation_distance(x, y, method = distance_method)
    num_landmarks <- ifelse(nrow(expression) > 1000, 500, nrow(expression))
    space <- SCORPIUS::reduce_dimensionality_landmarked(expression, dist_fun = dist_fun, ndim = ndim, num_landmarks = num_landmarks)
  }

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
      traj$path[-nrow(traj$path), , drop = FALSE] %>% set_colnames(., paste0("from_", colnames(.))),
      traj$path[-1, , drop = FALSE] %>% set_colnames(., paste0("to_", colnames(.)))
    )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotimes = traj$time
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
    mutate(time = prediction$pseudotimes[cell_id])

  seg_df <- prediction$dimred_trajectory_segments %>%
    as.data.frame

  g <- ggplot() +
    geom_path(aes(Comp1, Comp2), alpha = 0, data.frame(ranges %>% sweep(1, c(-.2, .2) * maxrange, "+"))) +
    geom_point(aes(Comp1, Comp2, colour = time), space_df) +
    geom_segment(aes(x = from_Comp1, xend = to_Comp1, y = from_Comp2, yend = to_Comp2), seg_df) +
    scale_colour_gradientn(colours = RColorBrewer::brewer.pal(3, "Dark2")) +
    labs(colour = "Pseudotime") +
    theme(legend.position = c(.92, .12))

  process_dynplot(g, prediction$id, expand = FALSE)
}
