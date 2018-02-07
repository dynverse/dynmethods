#' Description for Wishbone
#' @export
description_wishbone <- function() create_description(
  name = "Wishbone",
  short_name = "wishbone",
  package_loaded = c(),
  package_required = c("Wishbone"),
  par_set = makeParamSet(
    makeIntegerParam(id = "knn", lower = 2L, default = 15L, upper = 100L),
    makeIntegerParam(id = "n_diffusion_components", lower = 2L, default = 2L, upper = 20L),
    makeIntegerParam(id = "n_pca_components", lower = 2L, default = 15L, upper = 30L),
    makeLogicalParam(id = "branch", default = TRUE),
    makeIntegerParam(id = "k", lower = 2L, default = 15L, upper = 100L),
    makeIntegerParam(id = "num_waypoints", lower = 2L, default = 250L, upper = 500L),
    makeLogicalParam(id = "normalize", default = TRUE),
    makeNumericParam(id = "epsilon", lower = 0.1, default = 1, upper = 10)
  ),
  properties = c(),
  run_fun = run_wishbone,
  plot_fun = plot_wishbone
)

#' Description for Wanderlust
#' @export
description_wanderlust <- function() create_description(
  name = "Wanderlust",
  short_name = "wndrlust",
  package_loaded = c(),
  package_required = c("Wishbone"),
  par_set = makeParamSet(
    makeIntegerParam(id = "knn", lower = 2L, default = 15L, upper = 100L),
    makeIntegerParam(id = "n_diffusion_components", lower = 2L, default = 2L, upper = 20L),
    makeIntegerParam(id = "n_pca_components", lower = 2L, default = 15L, upper = 30L),
    makeLogicalParam(id = "branch", default = FALSE, tunable = FALSE),
    makeIntegerParam(id = "k", lower = 2L, default = 15L, upper = 100L),
    makeIntegerParam(id = "num_waypoints", lower = 2L, default = 250L, upper = 500L),
    makeLogicalParam(id = "normalize", default = TRUE),
    makeNumericParam(id = "epsilon", lower = 0.1, default = 1, upper = 10)
  ),
  properties = c(),
  run_fun = run_wishbone,
  plot_fun = plot_wishbone
)

run_wishbone <- function(
  counts,
  start_cells,
  marker_feature_ids = NULL,
  knn = 15,
  n_diffusion_components = 2,
  n_pca_components = 15,
  branch = TRUE,
  k = 15,
  num_waypoints = 250,
  normalize = TRUE,
  epsilon = 1
) {
  requireNamespace("Wishbone")

  if (num_waypoints > nrow(counts) / 2) {
    new_num <- ceiling(nrow(counts) / 2)
    warning("Reducing the number of waypoints to the number of cells (originally = ", num_waypoints, ", now = ", new_num, ")")
    num_waypoints <- new_num
  }

  start_cell_id <- sample(start_cells, 1)

  markers <-
    if (is.null(marker_feature_ids)) {
      "~"
    } else {
      marker_feature_ids
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

  # TIMING: after postproc
  tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network ,
    progressions = progressions,
    model = model
  ) %>% attach_timings_attribute(tl)
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
