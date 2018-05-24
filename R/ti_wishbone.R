abstract_wishbone_description <- function(method) {
  allow_branching <- method == "wishbone"
  name <- c("wishbone" = "Wishbone", "wndrlst" = "Wanderlust")[method] %>% setNames(NULL)

  create_ti_method(
    name = name,
    short_name = method,
    package_loaded = c(),
    package_required = c("Wishbone"),
    par_set = makeParamSet(
      makeIntegerParam(id = "knn", lower = 2L, default = 15L, upper = 100L),
      makeIntegerParam(id = "n_diffusion_components", lower = 2L, default = 2L, upper = 20L),
      makeIntegerParam(id = "n_pca_components", lower = 2L, default = 15L, upper = 30L),
      makeLogicalParam(id = "branch", default = allow_branching, tunable = allow_branching),
      makeIntegerParam(id = "k", lower = 2L, default = 15L, upper = 100L),
      makeIntegerParam(id = "num_waypoints", lower = 2L, default = 250L, upper = 500L),
      makeLogicalParam(id = "normalize", default = TRUE),
      makeNumericParam(id = "epsilon", lower = 0.1, default = 1, upper = 10),
      makeDiscreteParam(id = "method_name", values = c("wndrlst", "wishbone"), default = method, tunable = FALSE)
    ),
    run_fun = "run_wishbone",
    plot_fun = "plot_wishbone"
  )
}

#' Inferring trajectories with Wanderlust/Wishbone
#'
#' @inherit ti_identity description
#'
#' @param knn Number of nearest neighbours for diffusion map
#' @param n_diffusion_components Number of diffusion components
#' @param n_pca_components Number of pca components
#' @param branch Whether to find a branching (wishbone) or linear (wanderlust) trajectory
#' @param k Number of nearest neighbors for graph construction
#' @param num_waypoints Number of waypoints to sample
#' @param normalize Whether to normalize the data
#' @param epsilon Gaussian standard deviation for converting distances to affinities, for diffusion map
#'
#' @rdname wishbone
#'
#' @include wrapper_create_ti_method.R
#'
#' @export
ti_wishbone <- abstract_wishbone_description("wishbone")

#' @rdname wishbone
#' @export
ti_wndrlst <- abstract_wishbone_description("wndrlst")

run_wishbone <- function(
  counts,
  start_cells,
  n_end_states = NULL,
  marker_feature_ids = NULL,
  knn = 15,
  n_diffusion_components = 2,
  n_pca_components = 15,
  branch = TRUE,
  k = 15,
  num_waypoints = 250,
  normalize = TRUE,
  epsilon = 1,
  method_name = "wishbone"
) {
  requireNamespace("Wishbone")

  # reduce the number of waypoints if needed
  if (num_waypoints > nrow(counts) / 2) {
    new_num <- ceiling(nrow(counts) / 2)
    warning("Reducing the number of waypoints to the number of cells (originally = ", num_waypoints, ", now = ", new_num, ")")
    num_waypoints <- new_num
  }

  # sample one start cell if more are given
  start_cell_id <- sample(start_cells, 1)

  # transform marker features
  markers <-
    if (is.null(marker_feature_ids)) {
      "~"
    } else {
      marker_feature_ids
    }

  # if the number of end states are given, derive branch parameter from n_end_states
  if (method_name == "wishbone" && !is.null(n_end_states)) {
    branch <- n_end_states > 1
  }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # execute wishbone
  out <- Wishbone::Wishbone(
    counts = counts,
    start_cell_id = start_cell_id,
    knn = knn,
    n_diffusion_components = n_diffusion_components,
    n_pca_components = n_pca_components,
    markers = markers,
    branch = branch,
    k = k,
    num_waypoints = num_waypoints,
    normalize = normalize,
    epsilon = epsilon
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # retrieve model
  model <- out$branch_assignment %>%
    left_join(out$trajectory, by = "cell_id") %>%
    left_join(out$space, by = "cell_id")

  # create mapping between branch and from-to columns
  fromto <- if (branch) {
    data_frame(from = c("M1", "M2", "M2"), to = c("M2", "M3", "M4"), branch = c(1, 2, 3))
  } else {
    data_frame(from = "M1", to = "M2", branch = c(1, 2, 3))
  }

  # create network
  milestone_network <- model %>%
    group_by(branch) %>%
    summarise(length = diff(range(time))) %>%
    left_join(fromto, by = "branch") %>%
    mutate(directed = TRUE) %>%
    select(from, to, length, directed)

  # create progressions
  progressions <- model %>%
    left_join(fromto, by = "branch") %>%
    group_by(branch) %>%
    mutate(percentage = dynutils::scale_minmax(time)) %>%
    ungroup() %>%
    select(cell_id, from, to, percentage)

  # get the milestone names
  milestone_ids <- sort(unique(c(milestone_network$from, milestone_network$to)))

  # return output
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    model = model,
    divergence_regions = NULL
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom viridis scale_colour_viridis
plot_wishbone <- function(prediction) {
  g <- ggplot() +
    geom_point(aes(Comp1, Comp2, color = time), prediction$model %>% mutate_at(c("Comp1", "Comp2"), dynutils::scale_minmax)) +
    viridis::scale_colour_viridis() +
    labs(colour = "Trajectory") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}
