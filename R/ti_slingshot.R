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
  package_required = c("slingshot", "dyndimred"),
  parameters = list(
    ndim = list(
      type = "integer",
      default = 3L,
      upper = 20L,
      lower = 2L,
      description = "The number of dimensions"),
    nclus = list(
      type = "integer",
      default = 5L,
      upper = 40L,
      lower = 2L,
      description = "Number of clusters"),
    dimred = list(
      type = "discrete",
      default = "pca",
      values = c("pca", "mds", "tsne", "ica", "lle", "mds_sammon", "mds_isomds", "mds_smacof", "umap"),
      description = "A character vector specifying which dimensionality reduction method to use.\nSee \\link[dyndimred:dimred]{dyndimred::dimred} for the list of available dimensionality reduction methods."),

    shrink = list(
      type = "numeric",
      default = 1,
      upper = 1,
      lower = 0,
      description = "logical or numeric between 0 and 1, determines whether and how \nmuch to shrink branching lineages toward their average prior to the split."),
    reweight = list(
      type = "logical",
      default = TRUE,
      values = c("TRUE", "FALSE"),
      description = "logical, whether to allow cells shared between lineages to be\nreweighted during curve-fitting. If \\code{TRUE}, cells shared between \nlineages will be weighted by: distance to nearest curve / distance to\ncurve."),
    thresh = list(
      type = "numeric",
      default = -3,
      upper = 5,
      lower = -5,
      description = "numeric, determines the convergence criterion. Percent change in the total distance from cells to their projections along curves must be less than thresh. Default is 0.001, similar to principal.curve."),

    maxit = list(
      type = "integer",
      default = 10L,
      upper = 50L,
      lower = 0L,
      description = "numeric, maximum number of iterations, see principal.curve."),
    stretch = list(
      type = "numeric",
      default = 2,
      upper = 5,
      lower = 0,
      description = "numeric factor by which curves can be extrapolated beyond endpoints. Default is 2, see principal.curve."),
    smoother = list(
      type = "discrete",
      default = "smooth.spline",
      values = c("smooth.spline", "loess", "periodic.lowess"),
      description = "choice of scatter plot smoother. Same as principal.curve, but \"lowess\" option is replaced with \"loess\" for additional flexibility."),

    shrink.method = list(
      type = "discrete",
      default = "cosine",
      values = c("cosine", "tricube", "density"),
      description = "character denoting how to determine the appropriate\namount of shrinkage for a branching lineage. Accepted values are the same\nas for \\code{kernel} in \\code{\\link{density}} (default is \\code{\"cosine\"}),\nas well as \\code{\"tricube\"} and \\code{\"density\"}. See 'Details' for more.")
  ),
  run_fun = "dynmethods::run_slingshot",
  plot_fun = "dynmethods::plot_slingshot",
  apt_dependencies = c("libcgal-dev", "libglu1-mesa-dev", "libglu1-mesa-dev")
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
  space <- dyndimred::dimred(expr, method = dimred, ndim = ndim)

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

