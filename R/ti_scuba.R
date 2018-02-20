#' Description for SCUBA
#' @export
description_scuba <- function() create_description(
  name = "SCUBA",
  short_name = "scuba",
  package_loaded = c(),
  package_required = c("jsonlite", "readr", "SCUBA"),
  par_set = makeParamSet(
    makeLogicalParam(id = "rigorous_gap_stats", default = TRUE),
    makeIntegerParam(id = "N_dim", lower = 2L, upper = 3L, default = 2L), # limit to 3, limitation of sklearn tsne
    makeNumericParam(id = "low_gene_threshold", lower = 0, upper = 5, default = 1),
    makeNumericParam(id = "low_gene_fraction_max", lower = 0, upper = 1, default = 0.7),
    makeIntegerParam(id = "min_split", lower=1L, upper = 100L, default = 15L),
    makeNumericParam(id = "min_percentage_split", lower = 0, upper = 1, default = 0.25)
  ),
  properties = c(),
  run_fun = run_scuba,
  plot_fun = plot_scuba
)


run_scuba <- function(counts,
                      timecourse = NULL,
                      rigorous_gap_stats = TRUE,
                      N_dim = 2,
                      low_gene_threshold = 1,
                      low_gene_fraction_max = 0.7,
                      min_split = 15,
                      min_percentage_split = 0.25) {
  requireNamespace("SCUBA")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run scuba
  out <- SCUBA::SCUBA(
    counts = counts,
    rigorous_gap_stats = rigorous_gap_stats,
    N_dim = N_dim,
    low_gene_threshold = low_gene_threshold,
    low_gene_fraction_max = low_gene_fraction_max,
    min_split = min_split,
    min_percentage_split = min_percentage_split,
    timecourse = timecourse
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # get milestones
  unique_labs <- sort(unique(out$labels))
  milestone_fun <- function(x) paste0("milestone_", x)
  milestone_ids <- milestone_fun(unique_labs)

  # construct network
  milestone_network <- out$new_tree %>%
    select(from = `Parent cluster`, to = `Cluster ID`) %>%
    filter(to %in% unique_labs) %>%
    mutate(
      from = milestone_fun(from),
      to = milestone_fun(to),
      length = 1,
      directed = TRUE
    )

  # put cells on edges
  both_directions <- bind_rows(
    milestone_network %>% select(from, to) %>% mutate(label = from, percentage = 0),
    milestone_network %>% select(from, to) %>% mutate(label = to, percentage = 1)
  )
  progressions <- data_frame(
    cell_id = rownames(counts),
    label = milestone_fun(out$labels)
  ) %>%
    left_join(both_directions, by = "label") %>%
    group_by(cell_id) %>%
    arrange(desc(percentage)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-label)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_trajectory_to_wrapper(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    divergence_regions = NULL
  ) %>% add_timings_to_wrapper(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom grid arrow
plot_scuba <- function(prediction) {
  requireNamespace("igraph")

  #retrieve data
  cids <- prediction$cell_ids
  mids <- prediction$milestone_ids
  mnet <- prediction$milestone_network
  progs <- prediction$progressions

  # set colours
  mid_col <- setNames(seq_along(mids), mids)

  # perform dimred on network
  gr <- igraph::graph_from_data_frame(mnet, vertices = mids)
  lay <- igraph::layout_as_tree(gr, root = "milestone_0") %>%
    dynutils::scale_uniform()
  rownames(lay) <- mids
  colnames(lay) <- c("x", "y")

  # get clusters of cells
  labs <- progs %>% mutate(label = ifelse(percentage == 0, from, to)) %>% .$label

  # make plot of clusters
  mil_df <- data.frame(id = mids, lay)
  edge_df <- data.frame(
    mnet,
    from = lay[mnet$from,],
    to = lay[mnet$to,]
  )
  cel_df <- data.frame(
    row.names = NULL,
    id = cids,
    lab = labs,
    lay[labs,]
  )

  g <- ggplot() +
    geom_jitter(aes(y, x, colour = lab), cel_df, width = .03, height = .03) +
    geom_segment(aes(x = from.y, xend = to.y, y = from.x, yend = to.x), edge_df, arrow = grid::arrow()) +
    scale_colour_manual(values = mid_col) +
    theme(legend.position = "none")
  process_dynplot(g, prediction$id)
}
