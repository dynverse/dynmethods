#' Inferring trajectories with Control: shuffle
#'
#' This control method will return the milestone network of the provided
#' gold standard, but will shuffle the cell positions randomly.
#'
#' @param dummy_param This parameter does not do anything.
#'
#' @export
#' @importFrom dynplot plot_default
ti_shuffle <- create_ti_method(
  name = "Control: shuffle",
  short_name = "shuffle",
  package_loaded = c(),
  package_required = c(),
  trajectory_types = c("linear", "bifurcation", "convergence"),
  topology_inference = "free",
  type = "control_test",
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
    dummy_param = list(
      type = "numeric",
      default = 0.5,
      upper = 1,
      lower = 0,
      description = "Dummy parameter")
  ),
  run_fun = "dynmethods::run_shuffle",
  plot_fun = dynplot::plot_default
)

run_shuffle <- function(
  counts,
  task,
  dummy_param = .5
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # permute cell labels
  allcells <- rownames(counts)
  mapper <- setNames(sample(allcells), allcells)
  progressions <- task$progressions %>% mutate(
    cell_id = mapper[cell_id]
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = task$cell_ids
  ) %>% add_trajectory(
    milestone_ids = task$milestone_ids,
    milestone_network = task$milestone_network,
    progressions = progressions,
    divergence_regions = task$divergence_regions
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

