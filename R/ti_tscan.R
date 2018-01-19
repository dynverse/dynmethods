#' Description for TSCAN
#' @export
description_tscan <- function() create_description(
  name = "TSCAN",
  short_name = "TSCAN",
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
  properties = c(),
  run_fun = run_tscan,
  plot_fun = plot_tscan
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

  # project cells onto segments
  out <- project_cells_to_segments(
    cluster_network = cluster_network,
    cluster_space = cluster_space,
    sample_space = sample_space,
    sample_cluster = cds_clus$clusterid,
    num_segments_per_edge = 100,
    milestone_rename_fun = function(x) paste0("M", x)
  )

  # TIMING: after postproc
  tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(counts),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    space = out$space_df,
    centers = out$centers_df,
    edges = out$edge_df
  ) %>% attach_timings_attribute(tl)
}

plot_tscan <- function(prediction) {
  g <- ggplot() +
    geom_point(aes(PC1, PC2, colour = label, shape = label), prediction$space) +
    geom_text(aes(PC1, PC2, label = clus_id), prediction$centers, size = 6) +
    geom_segment(aes(x = from.PC1, xend = to.PC1, y = from.PC2, yend = to.PC2), prediction$edges) +
    theme(legend.position = "none")
  process_dynplot(g, prediction$id)
}

