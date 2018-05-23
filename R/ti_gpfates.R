#' Inferring trajectories with GPfates
#'
#' Arguments passed to this function will be used as default parameters for the method.
#'
#' @export
#'
#' @include wrapper_create_description.R
description_gpfates <- create_description(
  name = "GPfates",
  short_name = "gpfates",
  package_loaded = c(),
  package_required = c("GPfates"),
  par_set = makeParamSet(
    makeNumericParam(id = "log_expression_cutoff", lower = 0.5, upper = 5, default = 2),
    makeNumericParam(id = "min_cells_expression_cutoff", lower = 0, upper = 20, default = 2),
    makeIntegerParam(id = "ndims", lower = 1L, upper = 5L, default = 2L)
  ),
  run_fun = "run_gpfates",
  plot_fun = "plot_gpfates"
)

## TODO: give simulationtime as prior
#' @importFrom readr read_csv
#' @importFrom utils write.table
run_gpfates <- function(
  counts,
  n_end_states,
  ndims = 2,
  log_expression_cutoff = 2,
  min_cells_expression_cutoff = 2,
  num_cores = 1,
  verbose = FALSE
) {
  # documentation was not very detailed, so we had a hard time figuring out what the parameters were
  requireNamespace("GPfates")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # if the dataset is cyclic, pretend it isn't
  if (n_end_states == 0) {
    n_end_states <- 1
  }

  gp_out <- GPfates::GPfates(
    counts = counts,
    nfates = n_end_states,
    ndims = ndims,
    log_expression_cutoff = log_expression_cutoff,
    min_cells_expression_cutoff = min_cells_expression_cutoff,
    num_cores = num_cores,
    verbose = verbose
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  pseudotime <- gp_out$pseudotime %>% mutate(time = dynutils::scale_minmax(time))
  phi <- gp_out$phi
  dr <- gp_out$dr

  # get percentages of end milestones by multiplying the phi with the pseudotime
  progressions <- phi %>%
    gather(to, percentage, -cell_id) %>%
    left_join(pseudotime, by = "cell_id") %>%
    mutate(
      percentage = percentage*time,
      from = "M0"
    ) %>%
    select(cell_id, from, to, percentage)

  #  create milestone network
  milestone_network <- data_frame(
    from = "M0",
    to = paste0("M", seq_len(n_end_states)),
    length = 1,
    directed = TRUE
  )
  milestone_ids <- paste0("M", seq(0, n_end_states))

  # convert dimred and pseudotimes to right format
  dimred <- dr %>%
    as.data.frame() %>%
    magrittr::set_rownames(., .$cell_id) %>%
    select(-cell_id) %>%
    as.matrix()
  pseudotimes <- setNames(pseudotime$time, pseudotime$cell_id)

  divergence_regions <-
    data_frame(
      divergence_id = "divergence",
      milestone_id = paste0("M", seq(0, n_end_states)),
      is_start = milestone_id == "M0"
    )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    pseudotimes = pseudotimes,
    divergence_regions = divergence_regions
  ) %>% add_dimred(
    dimred = dimred
  ) %>% add_timings(
    tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_gpfates <- function(prediction, type = c("dimred", "assignment")) {
  type <- match.arg(type)

  sample_df <- prediction$dimred %>%
    as.data.frame %>%
    rownames_to_column("cell_id") %>%
    mutate(time = prediction$pseudotime[cell_id])

  switch(
    type,
    dimred = {
      max_trend <- prediction$milestone_percentages %>%
        group_by(cell_id) %>%
        arrange(desc(percentage)) %>%
        slice(1) %>%
        ungroup()

      plot_df <- sample_df %>%
        left_join(max_trend, by = "cell_id") %>%
        mutate(trend = gsub("M", "Trend ", milestone_id))

      g <- ggplot() +
        geom_point(aes(Comp1, Comp2, colour = trend), plot_df) +
        scale_color_brewer(palette = "Set2") +
        labs(colour = "Fitted trend") +
        theme(legend.position = c(.92, .12))

      process_dynplot(g, prediction$id)

      # TODO: Extract OGMP for plotting purposes. See dyneval #21
    },
    assignment = {
      progression_df <- sample_df %>%
        left_join(prediction$progressions, by = "cell_id") %>%
        mutate(trend = gsub("M", "Trend ", to))

      g <- ggplot() +
        geom_point(aes(time, percentage, colour = trend), progression_df) +
        facet_wrap(~trend, ncol = 1) +
        labs(x = "Pseudotime", y = "Trend assignment probability") +
        theme(legend.position = "none")

      process_dynplot(g, prediction$id)
    }
  )
}
