#' Description for SCOUP
#' @export
description_scoup <- function() create_description(
  name = "SCOUP",
  short_name = "SCOUP",
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
  properties = c(),
  run_fun = run_scoup,
  plot_fun = plot_scoup
)

#' @importFrom utils read.table write.table
#' @importFrom stats var
run_scoup <- function(
  expression,
  grouping_assignment,
  start_cells,
  n_end_states,
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

  start_cell <- sample(start_cells, 1)
  # figure out indices of starting population
  # from the grouping_assignment and the start_cell
  start_ix <- grouping_assignment %>%
    filter(cell_id %in% start_cell) %>%
    select(group_id) %>%
    left_join(grouping_assignment, by = "group_id") %>%
    .$cell_id

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run SP and SCOUP
  model <- SCOUP::run_SCOUP(
    expr = expression,
    start_ix = start_ix,
    ndim = ndim,
    nbranch = n_end_states,
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
  milestone_ids <- c("M0", paste0("M", seq_len(n_end_states)))

  # create milestone network
  milestone_network <- data_frame(
    from = milestone_ids[[1]],
    to = milestone_ids[-1],
    length = 1,
    directed = TRUE
  )

  # TIMING: after postproc
  tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    model = model
  ) %>% attach_timings_attribute(tl)
}

plot_scoup <- function(prediction, type = "dimred") {
  palette <- setNames(seq_along(prediction$milestone_ids), prediction$milestone_ids)

  # find most likely milestone of each cell
  celltypes <- prediction$milestone_percentages %>%
    group_by(cell_id) %>%
    arrange(desc(percentage)) %>%
    slice(1)

  # combine data
  space_df <- prediction$model$dimred %>%
    as.data.frame %>%
    rownames_to_column(var = "cell_id") %>%
    left_join(celltypes, by = "cell_id")

  # make plot
  g <- ggplot(space_df) +
    geom_point(aes(Comp1, Comp2, colour = milestone_id), shape = 1) +
    scale_colour_manual(values = palette) +
    labs(colour = "Milestone") +
    theme(legend.position = c(.92, .12))

  process_dynplot(g, id = prediction$id)
}
