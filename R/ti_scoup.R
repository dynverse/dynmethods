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
  package_required = c(),
  package_loaded = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "ndim", lower = 2L, default = 2L, upper = 20L),
    makeNumericParam(id = "max_ite1", lower = log(2), default = log(100), upper = log(5000), trafo = function(x) round(exp(x))), # should be 1000
    makeNumericParam(id = "max_ite2", lower = log(2), default = log(100), upper = log(50000), trafo = function(x) round(exp(x))), # should be 10000
    makeNumericParam(id = "alpha_min", lower = log(.001), default = log(.1), upper = log(10), trafo = exp),
    makeNumericParam(id = "alpha_max", lower = log(1), default = log(100), upper = log(10000), trafo = exp),
    makeNumericParam(id = "t_min", lower = log(.00001), default = log(.001), upper = log(1), trafo = exp),
    makeNumericParam(id = "t_max", lower = log(.1), default = log(2), upper = log(100), trafo = exp),
    makeNumericParam(id = "sigma_squared_min", lower = log(.001), default = log(.1), upper = log(10), trafo = exp),
    makeNumericParam(id = "thresh", lower = log(.01), default = log(.01), upper = log(10), trafo = exp)
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

  # create progressions
  milestone_percentages <- model$cpara %>%
    as.data.frame() %>%
    as_data_frame() %>%
    rownames_to_column("cell_id") %>%
    mutate(time = max(time) - time) %>%
    rename(M0 = time) %>%
    gather(milestone_id, percentage, -cell_id) %>%
    group_by(milestone_id) %>%
    mutate(percentage = percentage / max(percentage)) %>%
    ungroup() %>%
    group_by(cell_id) %>%
    mutate(percentage = percentage / sum(percentage)) %>%
    filter(percentage > 0 | milestone_id == "M0")

  # create milestone ids
  milestone_ids <- c("M0", paste0("M", seq_len(end_n)))

  # create milestone network
  milestone_network <- data_frame(
    from = milestone_ids[[1]],
    to = milestone_ids[-1],
    length = 1,
    directed = TRUE
  )

  # create divergence_regions
  divergence_regions <- data_frame(
    divergence_id = "divergence",
    milestone_id = milestone_ids,
    is_start = milestone_id == milestone_ids[[1]]
  )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression),
    cell_info = model$cpara %>% rownames_to_column("cell_id")
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    divergence_regions = divergence_regions
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
