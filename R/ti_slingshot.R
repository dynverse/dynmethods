#' Description for slingshot
#' @export
description_slingshot <- function() create_description(
  name = "slingshot",
  short_name = "slngsht",
  package_loaded = c(),
  package_required = c("slingshot"),
  par_set = makeParamSet(
    makeIntegerParam(id = "ndim", lower = 2L, upper = 20L, default = 3L),
    makeIntegerParam(id = "nclus", lower = 2L, upper = 40L, default = 5L),
    makeNumericParam(id = "shrink", lower = 0, upper = 1, default=1),
    makeLogicalParam(id = "reweight", default = TRUE),
    makeLogicalParam(id = "drop.multi", default = TRUE),
    makeNumericParam(id = "thresh", lower = -5, upper = 5, default = -3, trafo = function(x) 10^x),
    makeIntegerParam(id = "maxit", lower = 0L, upper = 50L, default = 10L),
    makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 2),
    makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "loess", "periodic.lowess")),
    makeDiscreteParam(id = "shrink.method", default = "cosine", values = c("cosine", "tricube", "density"))

  ),
  properties = c(),
  run_fun = run_slingshot,
  plot_fun = plot_slingshot
)

#' @importFrom stats prcomp kmeans
run_slingshot <- function(
  counts,
  start_cells = NULL,
  end_cells = NULL,
  ndim = 3,
  nclus = 5,
  dimred_name = "pca",
  shrink = 1,
  reweight = TRUE,
  drop.multi = TRUE,
  thresh = 0.001,
  maxit = 15,
  stretch = 2,
  smoother = "smooth.spline",
  shrink.method = "cosine"
) {
  requireNamespace("slingshot")

  start_cell <-
    if (!is.null(start_cells)) {
      sample(start_cells, 1)
    } else {
      NULL
    }

  # normalization & preprocessing
  # from the vignette of slingshot
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

  expr <- FQnorm(t(counts))

  # dimensionality reduction
  space <- stats::prcomp(t(log1p(expr)), scale. = FALSE)$x[,seq_len(ndim)]

  # clustering
  labels <- stats::kmeans(space, centers = nclus)$cluster

  # process prior data
  if(!is.null(start_cell)) {
    start.clus <- labels[[start_cell]]
  } else {
    start.clus <- NULL
  }
  if(!is.null(end_cells)) {
    end.clus <- unique(labels[end_cells])
  } else {
    end.clus <- NULL
  }

  # run slingshot
  sds <- slingshot::slingshot(
    reducedDim = space,
    clusterLabels = labels,
    start.clus = start.clus,
    end.clus = end.clus,
    shrink = shrink,
    reweight = reweight,
    drop.multi = drop.multi,
    thresh = thresh,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother,
    shrink.method = shrink.method
  )

  # adapted from plot-SlingshotDataSet
  # extract information on clusters
  lineages <- slingshot::lineages(sds)
  lineage_ctrl <- slingshot::lineageControl(sds)
  connectivity <- slingshot::connectivity(sds)
  clusterLabels <- slingshot::clusterLabels(sds)

  # calculate cluster centers
  centers <- t(sapply(rownames(connectivity), function(cli){
    colMeans(space[clusterLabels == cli,])
  }))

  # collect milestone network
  cluster_network <- lineages %>%
    map_df(~ data_frame(from = .[-length(.)], to = .[-1])) %>%
    unique() %>%
    mutate(
      length = lineage_ctrl$dist[cbind(from, to)],
      directed = TRUE # TODO: should be true
    )

  # project cells onto segments
  out <- project_cells_to_segments(
    cluster_network = cluster_network,
    cluster_space = centers,
    sample_space = sds@reducedDim,
    sample_cluster = clusterLabels,
    num_segments_per_edge = 100,
    milestone_rename_fun = function(x) paste0("M", x)
  )

  # collect curve data for visualisation purposes
  curves <- slingshot::curves(sds)
  curve_df <- names(curves) %>% map_df(function(id) {
    curve <- curves[[id]]
    data.frame(
      curve = id,
      curve$s,
      tag = curve$tag,
      lambda = curve$lambda,
      dist = curve$dist,
      w = curve$w,
      stringsAsFactors = FALSE
    )
  })

  # return output
  wrap_ti_prediction(
    ti_type = "dag?",
    id = "slingshot",
    cell_ids = rownames(counts),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    space = out$space_df,
    centers = out$centers_df,
    edge = out$edge_df,
    curve = curve_df
  )
}

#' @importFrom RColorBrewer brewer.pal
plot_slingshot <- function(prediction, type = c("lineage", "curve", "both")) {
  type <- match.arg(type)

  # reconstruct palette
  palt <- c(
    RColorBrewer::brewer.pal(9, "Set1")[-c(1,3,6)],
    RColorBrewer::brewer.pal(7, "Set2")[-2],
    RColorBrewer::brewer.pal(6, "Dark2")[-5],
    RColorBrewer::brewer.pal(8, "Set3")[-c(1,2)]
  )

  # apply to milestone ids
  cols <- setNames(palt[seq_along(prediction$milestone_ids)], prediction$milestone_ids)

  # create plots for curve if so requested
  if (type %in% c("curve", "both")) {
    gcurve <- geom_path(aes(PC1, PC2, group = curve), prediction$curve %>% arrange(curve, lambda))
  } else {
    gcurve <- NULL
  }

  # create plots for lineage if so requested
  if (type %in% c("lineage", "both")) {
    gcenter <- geom_point(aes(PC1, PC2), prediction$centers, size = 3)
    gsegment <- geom_segment(aes(x = from.PC1, xend = to.PC1, y = from.PC2, yend = to.PC2), prediction$edge)
  } else {
    gcenter <- NULL
    gsegment <- NULL
  }

  # return plot
  g <- ggplot() +
    gcurve +
    gsegment +
    geom_point(aes(PC1, PC2, colour = label), prediction$space) +
    gcenter +
    scale_colour_manual(values = cols) +
    labs(colour = "Milestone") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}

