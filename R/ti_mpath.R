#' Description for Mpath
#' @export
description_mpath <- function() create_description(
  name = "Mpath",
  short_name = "mpath",
  package_loaded = c("Mpath"),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "distMethod", default = "euclidean", values = c("pearson", "kendall", "spearman", "euclidean")),
    makeDiscreteParam(id = "method", default = "kmeans", values = c("kmeans", "diversity", "size", "diversity_size")),
    makeIntegerParam(id = "numcluster", lower = 3L, default = 11L, upper = 30L),
    makeLogicalParam(id = "numcluster_null", default = T),
    makeNumericParam(id = "diversity_cut", lower = .1, default = .6, upper = 1),
    makeNumericParam(id = "size_cut", lower = .01, default = .05, upper = 1)
  ),
  properties = c(),
  run_fun = run_mpath,
  plot_fun = plot_mpath
)

#' @importFrom utils write.table
#' @importFrom stats na.omit
#' @importFrom reshape2 melt
run_mpath <- function(
  counts,
  grouping_assignment,
  distMethod = "euclidean",
  method = "kmeans",
  numcluster = 11,
  numcluster_null = TRUE,
  diversity_cut = .6,
  size_cut = .05
) {
  requireNamespace("igraph")

  if (numcluster_null) {
    numcluster <- NULL
  }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # collect sample info
  sample_info <- grouping_assignment %>% rename(GroupID = group_id) %>% as.data.frame

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # designate landmarks
  landmark_cluster <- Mpath::landmark_designation(
    rpkmFile = t(counts),
    baseName = NULL,
    sampleFile = sample_info,
    distMethod = distMethod,
    method = method,
    numcluster = numcluster,
    diversity_cut = diversity_cut,
    size_cut = size_cut,
    saveRes = FALSE
  ) %>%
    mutate_if(is.factor, as.character)

  milestone_ids <- unique(landmark_cluster$landmark_cluster)

  # catch situation where mpath only detects 1 landmark
  if (length(milestone_ids) == 1) {
    stop("Mpath only detected one landmark")
  }

  # build network
  network <- Mpath::build_network(
    exprs = t(counts),
    baseName = NULL,
    landmark_cluster = landmark_cluster,
    distMethod = distMethod,
    writeRes = FALSE
  )

  # trim network
  trimmed_network <- Mpath::trim_net(
    nb12 = network,
    writeRes = FALSE
  )

  # create final milestone network
  class(trimmed_network) <- NULL
  milestone_network <- trimmed_network %>%
    reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
    mutate_if(is.factor, as.character) %>%
    filter(length > 0, from < to) %>%
    mutate(directed = FALSE)

  grouping <-
    with(landmark_cluster, setNames(landmark_cluster, cell))

  wrap_prediction_model(
    cell_ids = rownames(counts),
    grouping_assignment = grouping_assignment
  ) %>% add_grouping_to_wrapper(
    group_ids = milestone_ids,
    grouping = grouping
  ) %>% add_cluster_graph_to_wrapper(
    milestone_network = milestone_network
  ) %>% add_timings_to_wrapper(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_mpath <- function(prediction) {
  requireNamespace("igraph")
  requireNamespace("ggforce")
  requireNamespace("RColorBrewer")

  # milestone net as igraph in order to perform dimred
  edges <- prediction$milestone_network %>% filter(to != "FILTERED_CELLS")
  gr <- igraph::graph_from_data_frame(
    edges,
    directed = FALSE,
    vertices = prediction$milestone_ids
  )
  lay <- gr %>%
    igraph::layout_with_kk() %>%
    dynutils::scale_quantile(0)
  dimnames(lay) <- list(
    igraph::V(gr)$name,
    c("X", "Y")
  )
  lay_df <- lay %>% as.data.frame %>% rownames_to_column("milestone_id")

  # collect info on cells
  cell_ids <- prediction$cell_ids
  labels <- prediction$grouping_assignment %>%
    slice(match(cell_ids, cell_id)) %>%
    .$group_id
  clustering <- prediction$progressions %>%
    slice(match(cell_ids, cell_id)) %>%
    {with(., ifelse(percentage == 0, from, to))}

  # generate pie df with positioning
  pie_df <- data_frame(cell_id = cell_ids, label = labels, milestone_id = clustering) %>%
    group_by(milestone_id, label) %>%
    summarise(n = n()) %>%
    mutate(
      value = n / sum(n) * 2 * pi,
      end = cumsum(value),
      start = end - value
    ) %>%
    ungroup() %>%
    left_join(lay_df, by = "milestone_id")

  # generate edge df with positioning
  edges_df <- edges %>%
    left_join(lay_df %>% select(from = milestone_id, from.x = X, from.y = Y), by = "from") %>%
    left_join(lay_df %>% select(to = milestone_id, to.x = X, to.y = Y), by = "to")

  # Determine a colour scheme
  ann_groups <- unique(labels) %>% sort
  ann_cols <- setNames(RColorBrewer::brewer.pal(length(ann_groups), "Set2"), ann_groups)

  # Make a line plot
  g <- ggplot() +
    geom_segment(aes(x = from.x, xend = to.x, y = from.y, yend = to.y), edges_df) +
    ggforce::geom_arc_bar(aes(x0 = X, y0 = Y, r0 = 0, r = .075,
                              start = start, end = end, fill = label, group = milestone_id), data = pie_df) +
    geom_text(aes(X, Y, label = milestone_id), lay_df) +
    scale_fill_manual(values = ann_cols) +
    theme(legend.position = c(.92, .12))

  process_dynplot(g, prediction$id)
}


