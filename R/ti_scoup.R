#' Inferring trajectories with SCOUP
#'
#' @inherit ti_angle description
#'
#' @param ndim Number of pca dimensions
#' @param max_ite1 Upper bound of EM iteration (without pseudo-time optimization). The detailed explanation is described in the supplementary text. (default is 1,000)
#' @param max_ite2 Upper bound of EM iteration (including pseudo-time optimization) (default is 1,000).
#' @param alpha_min Lower bound of alpha (default is 0.1)
#' @param alpha_max Upper bound of alpha (default is 100)
#' @param t_min Lower bound of pseudo-time (default is 0.001)
#' @param t_max Upper bound of pseudo-time (default is 2.0)
#' @param sigma_squared_min Lower bound of sigma squared (default is 0.1)
#' @param thresh Threshold
#'
#' @export
ti_scoup <- create_ti_method(
  name = "SCOUP",
  short_name = "scoup",
  package_required = c("SCOUP"),
  package_loaded = c(),
  doi = "10.1186/s12859-016-1109-3",
  trajectory_types = "linear",
  topology_inference = "parameter",
  type = "algorithm",
  license = "MIT",
  authors = list(
    list(
      given = "Hirotaka",
      family = "Matsumoto",
      email = "hirotaka.matsumoto@riken.jp",
      github = "hmatsu1226"
    )
  ),
  publication_date = "2016-06-08",
  code_url = "https://github.com/gcyuan/PySCUBA",
  parameters = list(
    ndim = list(
      type = "integer",
      default = 2L,
      upper = 20L,
      lower = 2L,
      description = "Number of pca dimensions"
    ),
    max_ite1 = list(
      type = "numeric",
      default = 100,
      upper = 5000,
      lower = 2,
      description = "Upper bound of EM iteration (without pseudo-time optimization). The detailed explanation is described in the supplementary text. (default is 1,000)"
    ),
   max_ite2 = list(
      type = "numeric",
      default = 100,
      upper = 500000,
      lower = 2,

      description = "Upper bound of EM iteration (including pseudo-time optimization) (default is 1,000)."
    ),
    alpha_min = list(
      type = "numeric",
      default = 0.1,
      upper = 10,
      lower = 0.001,
      description = "Lower bound of alpha (default is 0.1)"
    ),
    alpha_max = list(
      type = "numeric",
      default = 100,
      upper = 10000,
      lower = 1,
      description = "Upper bound of alpha (default is 100)"
    ),
    t_min = list(
      type = "numeric",
      default = 0.001,
      upper = 1,
      lower = 0.00001,
      description = "Lower bound of pseudo-time (default is 0.001)"
    ),
    t_max = list(
      type = "numeric",
      default = 2,
      upper = 100,
      lower = 0.1,
      description = "Upper bound of pseudo-time (default is 2.0)"
    ),
    sigma_squared_min = list(
      type = "numeric",
      default = 0.1,
      upper = 10,
      lower = 0.001,
      description = "Lower bound of sigma squared (default is 0.1)"
    ),
    thresh = list(
      type = "numeric",
      default = 0.01,
      upper = 10,
      lower = 0.01,
      description = "Threshold"
    )
  ),
  run_fun = "dynmethods::run_scoup",
  plot_fun = "dynmethods::plot_scoup"
)

#' @importFrom utils read.table write.table
#' @importFrom stats var
run_scoup <- function(
  expression,
  groups_id,
  start_id,
  end_n,
  ndim = 2,
  max_ite1 = 100,
  max_ite2 = 100,
  alpha_min = .1,
  alpha_max = 100,
  t_min = .001,
  t_max = 2,
  sigma_squared_min = .1,
  thresh = .01,
  verbose = FALSE
) {
  requireNamespace("SCOUP")

  # if the dataset is cyclic, pretend it isn't
  if (end_n == 0) {
    end_n <- 1
  }

  start_cell <- sample(start_id, 1)
  # figure out indices of starting population
  # from the groups_id and the start_cell
  start_ix <- groups_id %>%
    filter(cell_id %in% start_cell) %>%
    select(group_id) %>%
    left_join(groups_id, by = "group_id") %>%
    .$cell_id

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run SP and SCOUP
  model <- SCOUP::run_SCOUP(
    expr = expression,
    start_ix = start_ix,
    ndim = ndim,
    nbranch = end_n,
    max_ite1 = max_ite1,
    max_ite2 = max_ite2,
    alpha_min = alpha_min,
    alpha_max = alpha_max,
    t_min = t_min,
    t_max = t_max,
    sigma_squared_min = sigma_squared_min,
    thresh = thresh,
    verbose = verbose
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  if (any(is.na(model$ll))) {
    stop("SCOUP returned NaNs", call. = FALSE)
  }

  pseudotime <- model$cpara %>% {set_names(.$time, rownames(.))}
  esp <- model$cpara %>% select(-time) %>% rownames_to_column("cell_id")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression),
    cell_info = model$cpara %>% tibble::rownames_to_column("cell_id")
  ) %>% add_end_state_probabilities(
    end_state_probabilities = esp,
    pseudotime = pseudotime,
    do_scale_minmax = TRUE
  ) %>% add_dimred(
    dimred = model$dimred %>% as.matrix
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_scoup <- function(prediction, type = "dimred") {
  palette <- setNames(seq_along(prediction$milestone_ids), prediction$milestone_ids)

  # find most likely milestone of each cell
  celltypes <- prediction$milestone_percentages %>%
    group_by(cell_id) %>%
    arrange(desc(percentage)) %>%
    slice(1)

  # combine data
  space_df <- prediction$dimred %>%
    as.data.frame %>%
    rownames_to_column(var = "cell_id") %>%
    left_join(celltypes, by = "cell_id")

  # make plot
  g <- ggplot(space_df) +
    geom_point(aes(comp_1, comp_2, colour = milestone_id), shape = 1) +
    scale_colour_manual(values = palette) +
    labs(colour = "Milestone") +
    theme(legend.position = c(.92, .12))

  process_dynplot(g, id = prediction$id)
}
