#' Inferring trajectories with Slingshot
#'
#' @inherit ti_angle description
#'
#' @param nclus Number of clusters
#' @param dimred A character vector specifying which dimensionality reduction method to use.
#'   See [dyndimred::dimred] for the list of available dimensionality reduction methods.
#' @param thresh numeric, determines the convergence criterion. Percent change in the total distance from cells to their projections along curves must be less than thresh. Default is 0.001, similar to principal.curve.
#' @param maxit numeric, maximum number of iterations, see principal.curve.
#' @param stretch numeric factor by which curves can be extrapolated beyond endpoints. Default is 2, see principal.curve.
#' @param smoother choice of scatter plot smoother. Same as principal.curve, but "lowess" option is replaced with "loess" for additional flexibility.
#' @inheritParams dyndimred::dimred
#' @inheritParams slingshot::slingshot
#'
#' @seealso [slingshot::slingshot()]
#'
#' @export
ti_slingshot <- create_ti_method(
  name = "Slingshot",
  short_name = "slingshot",
  package_loaded = c(),
  package_required = c("slingshot"),
  par_set = makeParamSet(
    makeIntegerParam(id = "ndim", lower = 2L, upper = 20L, default = 3L),
    makeIntegerParam(id = "nclus", lower = 2L, upper = 40L, default = 5L),
    makeDiscreteParam(id = "dimred", default = "pca", values = names(dyndimred::list_dimred_methods())),
    makeNumericParam(id = "shrink", lower = 0, upper = 1, default=1),
    makeLogicalParam(id = "reweight", default = TRUE),
    makeLogicalParam(id = "drop.multi", default = TRUE),
    makeNumericParam(id = "thresh", lower = -5, upper = 5, default = -3, trafo = function(x) 10^x),
    makeIntegerParam(id = "maxit", lower = 0L, upper = 50L, default = 10L),
    makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 2),
    makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "loess", "periodic.lowess")),
    makeDiscreteParam(id = "shrink.method", default = "cosine", values = c("cosine", "tricube", "density"))

  ),
  run_fun = "dynmethods::run_slingshot",
  plot_fun = "dynmethods::plot_slingshot"
)

run_slingshot <- function(
  counts,
  start_id = NULL,
  end_id = NULL,
  ndim = 3,
  nclus = 5,
  dimred = "pca",
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
    if (!is.null(start_id)) {
      sample(start_id, 1)
    } else {
      NULL
    }

  # normalization & preprocessing
  # from the vignette of slingshot
  FQnorm <- function(counts){
    rk <- apply(counts, 2, rank, ties.method = "min")
    counts.sort <- apply(counts, 2, sort)
    refdist <- apply(counts.sort, 1, median)
    norm <- apply(rk, 2, function(r) refdist[r])
    rownames(norm) <- rownames(counts)
    return(norm)
  }

  expr <- t(log1p(FQnorm(t(counts))))

  # dimensionality reduction
  space <- dimred(expr, method = dimred, ndim = ndim)

  # clustering
  labels <- stats::kmeans(space, centers = nclus)$cluster

  # process prior data
  if(!is.null(start_cell)) {
    start.clus <- labels[[start_cell]]
  } else {
    start.clus <- NULL
  }
  if(!is.null(end_id)) {
    end.clus <- unique(labels[end_id])
  } else {
    end.clus <- NULL
  }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run slingshot
  sds <- slingshot::slingshot(
    space,
    labels,
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

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # adapted from plot-SlingshotDataSet
  # extract information on clusters
  lineages <- slingshot::slingLineages(sds)
  lineage_ctrl <- slingshot::slingParams(sds)
  connectivity <- slingshot::slingAdjacency(sds)
  clusterLabels <- slingshot::clusterLabels(sds) %>% setNames(rownames(counts))

  # calculate cluster centers
  centers <- t(sapply(rownames(connectivity), function(cli){
    colMeans(space[clusterLabels[, cli] == 1,,drop=T])
  }))

  # collect milestone network
  cluster_network <- lineages %>%
    map_df(~ data_frame(from = .[-length(.)], to = .[-1])) %>%
    unique() %>%
    mutate(
      length = lineage_ctrl$dist[cbind(from, to)],
      directed = TRUE # TODO: should be true
    )

  # collect curve data for visualisation purposes
  curves <- slingshot::slingCurves(sds)
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
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_dimred_projection(
    milestone_ids = rownames(centers),
    milestone_network = cluster_network,
    dimred_milestones = centers,
    dimred = sds@reducedDim,
    milestone_assignment_cells = clusterLabels,
    num_segments_per_edge = 100,
    curve = curve_df
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_slingshot <- function(prediction, type = c("curve", "lineage", "both")) {
  requireNamespace("RColorBrewer")

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
    gcurve <- geom_path(aes(comp_1, comp_2, group = curve), prediction$curve %>% arrange(curve, lambda))
  } else {
    gcurve <- NULL
  }

  # create plots for lineage if so requested
  if (type %in% c("lineage", "both")) {
    gcenter <- geom_point(aes(comp_1, comp_2), prediction$dimred_milestones %>% as.data.frame, size = 3)
    gsegment <- geom_segment(aes(x = from_comp_1, xend = to_comp_1, y = from_comp_2, yend = to_comp_2), prediction$dimred_trajectory_segments %>% as.data.frame())
  } else {
    gcenter <- NULL
    gsegment <- NULL
  }

  space <- prediction$dimred %>%
    as.data.frame %>%
    rownames_to_column("cell_id") %>%
    mutate(label = prediction$milestone_assignment_cells[cell_id])

  # return plot
  g <- ggplot() +
    geom_point(aes(comp_1, comp_2, colour = label), space) +
    gcurve +
    gsegment +
    gcenter +
    scale_colour_manual(values = cols) +
    labs(colour = "Milestone") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}

