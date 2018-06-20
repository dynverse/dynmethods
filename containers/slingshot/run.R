library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(slingshot)
library(dyndimred)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  counts,
  start_id = NULL,
  end_id = NULL,
  ndim = 3,
  nclus = 5,
  dimred = "pca",
  shrink = 1,
  reweight = TRUE,
  thresh = -3,
  maxit = 10,
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

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')