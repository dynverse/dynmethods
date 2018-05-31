#' Inferring trajectories with TSCAN
#'
#' @inherit ti_angle description
#'
#' @inheritParams TSCAN::preprocess
#' @inheritParams TSCAN::exprmclust
#' @param clusternum_lower An integer specifying the minimal possible cluster numbers. The best cluster number will be picked using BIC.
#' @param clusternum_upper An integer specifying the maximal possible cluster numbers. The best cluster number will be picked using BIC.
#'
#' @export
ti_tscan <- create_ti_method(
  name = "TSCAN",
  short_name = "tscan",
  package_loaded = c(),
  package_required = c("TSCAN", "igraph"),
  par_set = makeParamSet(
    makeNumericParam(id = "minexpr_percent", lower=0, upper=1, default=0),
    makeNumericParam(id = "minexpr_value", lower=0, upper=10, default=0),
    makeNumericParam(id = "cvcutoff", lower=0, upper=5, default=0),
    makeIntegerParam(id = "clusternum_lower", lower = 2L, upper = 20L, default = 2L),
    makeIntegerParam(id = "clusternum_upper", lower = 2L, upper = 20L, default = 9L),
    makeDiscreteParam(id = "modelNames", default = "VVV", values = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV")),
    forbidden = quote(clusternum_lower > clusternum_upper)
  ),
  run_fun = "dynmethods::run_tscan",
  plot_fun = "dynmethods::plot_tscan"
)

run_tscan <- function(
  counts,
  minexpr_percent = 0,
  minexpr_value = 0,
  cvcutoff = 0,
  clusternum_lower = 2,
  clusternum_upper = 9,
  modelNames = "VVV"
) {
  requireNamespace("TSCAN")
  requireNamespace("igraph")

  # process clusternum
  clusternum <- seq(clusternum_lower, clusternum_upper, 1)

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # preprocess counts
  cds_prep <- TSCAN::preprocess(
    t(as.matrix(counts)),
    takelog = TRUE,
    logbase = 2,
    pseudocount = 1,
    clusternum = NULL,
    minexpr_value = minexpr_value,
    minexpr_percent = minexpr_percent,
    cvcutoff = cvcutoff
  )

  # cluster the data
  cds_clus <- TSCAN::exprmclust(
    cds_prep,
    clusternum = clusternum,
    modelNames = modelNames,
    reduce = TRUE
  )

  # order the cells
  cds_order <- TSCAN::TSCANorder(cds_clus)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # process output
  cluster_network <- cds_clus$MSTtree %>%
    igraph::as_data_frame() %>%
    rename(length = weight) %>%
    mutate(directed = FALSE)
  sample_space <- cds_clus$pcareduceres
  cluster_space <- cds_clus$clucenter
  rownames(cluster_space) <- as.character(seq_len(nrow(cluster_space)))
  colnames(cluster_space) <- colnames(sample_space)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_dimred_projection(
    milestone_ids = rownames(cluster_space),
    milestone_network = cluster_network,
    dimred_milestones = cluster_space,
    dimred = sample_space,
    milestone_assignment_cells = cds_clus$clusterid,
    num_segments_per_edge = 100
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_tscan <- function(prediction) {
  requireNamespace("ggrepel")

  space <-
    prediction$dimred %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    mutate(label = prediction$milestone_assignment_cells[cell_id])

  g <- ggplot() +
    geom_point(aes(PC1, PC2, colour = label, shape = label), space) +
    ggrepel::geom_text_repel(aes(PC1, PC2, label = clus_id), prediction$dimred_milestones %>% as.data.frame %>% rownames_to_column("clus_id") %>% mutate(clus_id = paste0("M", row_number())), size = 6, min.segment.length = Inf) +
    geom_segment(aes(x = from_PC1, xend = to_PC1, y = from_PC2, yend = to_PC2), prediction$dimred_trajectory_segments %>% as.data.frame()) +
    theme(legend.position = "none")
  process_dynplot(g, prediction$id)
}

