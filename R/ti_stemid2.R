#' Inferring trajectories with StemID2
#'
#' @param clustnr maximum number of clusters for the computation of the gap statistic or derivation of the cluster number by the saturation criterion. Default is 30. If more major cell types are expected a higher number should be chosen.
#' @param bootnr number of booststrapping runs for clusterboot. Default is 50
#' @param metric the input data are transformed to a distance object. Distances can be computed based on different metrics. Possible values are "pearson", "spearman", "logpearson", "euclidean", "kendall", "maximum", "manhattan", "canberra", "binary" or "minkowski". Default is "pearson". In case of the correlation based methods, the distance is computed as 1 – correlation. K-medoids clustering is performed on this distance object.
#' @param num_cluster_method the type of clustering method, can be sat, gap or manual
#' @param SE.method the clustering routine calls a modified version of the maxSE function from the cluster package to determine the first local maximum of the gap statistic. By default, we use the method "Tibs2001SEmax" for calling the first local maximum (see specification of maxSE). This method requires that the maximum exceeds the values of its neighbors by a fraction of their standard deviation. This fraction is defined by the parameter SE.factor. All methods defined for the original maxSE function can also be used.
#' @param SE.factor fraction of the standard deviation by which the local maximum is required to differ from the neighboring points it is compared to. Default is 0.25.
#' @param B.gap number of bootstrap runs for the calculation of the gap statistic. Default is 50.
#' @param cln the number of clusters for k-medoids clustering. Default is 0. In this case, the cluster number is determined based on the gap statistic6 and do.gap has to be TRUE.
#' @param FUNcluster the clustering method applied. One of the following methods can be selected: kmedoids, kmeans, hclust. RaceID3 is designed for k-medoids clustering and therefore it is recommended to use only the kmedoids method. Default is kmedoids.
#' @param dimred_method the dimensionality reduction method, can be tsne, sammon or tsne_initcmd
#' @param outminc expression cutoff for the identification of outlier genes is defined. Default is 5.
#' @param outlg minimal number of outlier genes required to identify a cell as an outlier. Default is 2.
#' @param probthr defines the probability threshold for outlier calling. If the probability of observing a given expression level for a gene in a cell is lower than this cutoff (based on the negative binomial distribution for the calibrated noise model), the cell is considered an outlier for this gene. Default is 10-3.
#' @param thr_lower lower probability for which the number of outliers is computed in order to plot the dependence of the number of outliers on the probability threshold
#' @param thr_upper upper probability for which the number of outliers is computed in order to plot the dependence of the number of outliers on the probability threshold
#' @param outdistquant outlier cells are merged to outlier clusters if their similarity exceeds the outdistquant-quantile of the similarity distribution for all pairs of cells that are together in one of the original clusters. Default is 0.95.
#' @param nmode Boolean argument. If nmode is set to TRUE the assignment to inter-cluster links for each cell is not done based on the longest projection, but based on identifying the cluster (other than the cluster the cell belongs to) that contains the nearest neighbor of the cell, i. e. the cell with the most similar transcriptome. The coordinate on the assigned link is still derived based on the projection. Default is FALSE.
#' @param pdishuf positive integer. This is the number of randomizations to be performed. As a rule of thumb this number should be at least one order of magnitude larger than the desired p-value on the significance of the number of cells on a connection. Default is 2000.
#' @param pthr positive number. This number corresponds to the p-value threshold, which is used to determine, whether the magnitude of an observed trajectory is significantly larger than observed for the randomized background distribution. This criterion is not used to infer significance of a link, but shown in a graphical representation of the tree
#' @param pethr positive number. This number corresponds to the p-value threshold, which isused to determine for each link if it is populated by a number of cellssignificantly larger than expected for the randomized background distribution. This p-value threshold determines, which connections are considered validdifferentiation trajectories in the derived lineage tree.
#'
#' @inherit ti_identity description
#'
#' @export
#'
#' @include wrapper_create_ti_method.R
ti_stemid2 <- create_ti_method(
  name = "StemID2",
  short_name = "stemid2",
  package_loaded = c(),
  package_required = c("StemID2"),
  par_set = makeParamSet(
    makeIntegerParam(id = "clustnr", lower = 10L, default = 30L, upper = 100L),
    makeIntegerParam(id = "bootnr", lower = 20L, default = 50L, upper = 100L),
    makeDiscreteParam(id = "metric", default = "pearson", values = c("pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")),
    makeDiscreteParam(id = "num_cluster_method", default = "sat", values = c("sat", "gap", "manual")),
    makeDiscreteParam(id = "SE.method", default = "Tibs2001SEmax", values = c("firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax")),
    makeNumericParam(id = "SE.factor", default = .25, lower = 0, upper = 1),
    makeIntegerParam(id = "B.gap", lower = 20L, default = 50L, upper = 100L),
    makeIntegerParam(id = "cln", lower = 20L, default = 30L, upper = 100L),
    makeDiscreteParam(id = "FUNcluster", default = "kmedoids", values = c("kmedoids", "kmeans", "hclust")),
    makeDiscreteParam(id = "dimred_method", default = "tsne", values = c("tsne", "sammon", "tsne_initcmd")),
    makeIntegerParam(id = "outminc", lower = 0L, default = 0L, upper = 100L), # default should be 5, but stemid otherwise frequently produces errors
    makeIntegerParam(id = "outlg", lower = 0L, default = 2L, upper = 100L),
    makeNumericParam(id = "probthr", lower = -10, default = -10, upper = -1, trafo = function(x) 10^x), # lowered to almost zero, was 10^-3
    makeNumericParam(id = "thr_lower", lower = -100, default = -40, upper = -1),
    makeNumericParam(id = "thr_upper", lower = -100, default = -1, upper = -1),
    makeNumericParam(id = "outdistquant", lower = 0, default = .95, upper = 1),
    makeLogicalParam(id = "nmode", default = FALSE),
    makeNumericParam(id = "pdishuf", lower = log(100), default = log(500), upper = log(10000), trafo = function(x) round(exp(x))), # orig 2000
    makeNumericParam(id = "pthr", lower = -4, default = -2, upper = 0, trafo = function(x) 10^x),
    makeNumericParam(id = "pethr", lower = -4, default = -2, upper = 0, trafo = function(x) 10^x),
    forbidden = quote(thr_lower > thr_upper)
  ),
  run_fun = "run_stemid2",
  plot_fun = "plot_stemid2"
)

run_stemid2 <- function(
  expression,
  clustnr = 30,
  bootnr = 50,
  metric = "pearson",
  num_cluster_method = "sat",
  SE.method = "Tibs2001SEmax",
  SE.factor = .25,
  B.gap = 50,
  cln = 30L,
  FUNcluster = "kmedoids",
  dimred_method = "tsne", # tsne, sammon, tsne_initcmd
  outminc = 0,
  outlg = 2,
  probthr = 1e-3,
  thr_lower = -40,
  thr_upper = -1,
  outdistquant = .95,
  nmode = FALSE,
  pdishuf = 2000,
  pthr = .01,
  pethr = .01
) {
  requireNamespace("StemID2")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # initialize SCseq object with transcript expression
  sc <- StemID2::SCseq(data.frame(t(expression), check.names = FALSE))

  # filtering of expression data
  sc <- sc %>% StemID2::filterdata(
    mintotal = 1,
    minexpr = 0,
    minnumber = 0,
    maxexpr = Inf,
    downsample = FALSE,
    sfn = FALSE, # newly added
    hkn = FALSE, # newly added
    dsn = 1,
    CGenes = NULL, # newly added
    FGenes = NULL, # newly added
    ccor = .4 # newly added
  )

  # k-medoids clustering
  do_gap <- num_cluster_method == "gap"
  do_sat <- num_cluster_method == "sat"

  sc <- tryCatch({
    sc %>% StemID2::clustexp(
      clustnr = clustnr,
      bootnr = bootnr,
      metric = metric,
      do.gap = do_gap,
      sat = do_sat,
      SE.method = SE.method,
      SE.factor = SE.factor,
      B.gap = B.gap,
      cln = cln,
      FUNcluster = FUNcluster,
      FSelect = TRUE
    )
  }, error = function(e) {
    sc %>% StemID2::clustexp(
      clustnr = clustnr,
      bootnr = bootnr,
      metric = metric,
      do.gap = do_gap,
      sat = do_sat,
      SE.method = SE.method,
      SE.factor = SE.factor,
      B.gap = B.gap,
      cln = cln,
      FUNcluster = FUNcluster,
      FSelect = FALSE # turn off feature filtering if clusterexp errors because of it
    )
  })

  # compute t-SNE map
  sammonmap <- dimred_method == "sammon"
  initial_cmd <- dimred_method == "tsne_initcmd"
  sc <- sc %>% StemID2::comptsne(
    sammonmap = sammonmap,
    initial_cmd = initial_cmd,
    fast = TRUE, # newly added
    perplexity = 30 # newly added
  )

  # detect outliers and redefine clusters
  sc <- sc %>% StemID2::findoutliers(
    outminc = 5,
    outlg = outlg,
    probthr = probthr,
    thr = 2^(thr_lower:thr_upper),
    outdistquant = outdistquant
  )

  sc <- sc %>% StemID2::rfcorrect(
    final = TRUE, # newly added
    nbfactor = 5 # newly added
  )

  # initialization
  ltr <- StemID2::Ltree(sc)

  # computation of the entropy
  ltr <- ltr %>% StemID2::compentropy()

  # computation of the projections for all cells
  ltr <- ltr %>% StemID2::projcells(
    cthr = 0, # default = 2
    nmode = nmode
  )

  # computation of the projections for all cells after randomization
  ltr <- ltr %>% StemID2::projback(
    pdishuf = pdishuf,
    nmode = nmode,
    fast = FALSE # newly added
  )

  # assembly of the lineage tree
  ltr <- ltr %>% StemID2::lineagetree(
    pthr = pthr,
    nmode = nmode,
    fast = FALSE # newly added
  )

  # compute a spanning tree
  ltr <- ltr %>% StemID2::compspantree()

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # get network info
  cluster_network <- data_frame(
    from = as.character(ltr@ldata$m[-1]),
    to = as.character(ltr@trl$trl$kid),
    length = ltr@trl$dc[cbind(from, to)],
    directed = FALSE
  )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_dimred_projection(
    milestone_ids = rownames(ltr@ldata$cnl %>% as.matrix),
    milestone_network = cluster_network,
    dimred_milestones = ltr@ldata$cnl %>% as.matrix,
    dimred = ltr@ltcoord,
    milestone_assignment_cells = as.character(ltr@ldata$lp) %>% setNames(rownames(expression)),
    num_segments_per_edge = 100,
    col_ann = ltr@sc@fcol
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_stemid2 <- function(prediction) {
  col_ann <- setNames(ltr@sc@fcol, prediction$milestone_ids)

  space <- prediction$dimred %>%
    as.data.frame %>%
    rownames_to_column("cell_id") %>%
    mutate(label = prediction$milestone_assignment_cells[cell_id])

  space_clus <- prediction$dimred_milestones %>%
    as.data.frame %>%
    rownames_to_column("clus_id")

  g <- ggplot() +
    geom_point(aes(V1, V2), space, size = 2, colour = "darkgray", na.rm = TRUE) +
    geom_text(aes(V1, V2, label = label, colour = label), space, na.rm = TRUE) +
    geom_text(aes(V1, V2, label = clus_id), space_clus, size = 8) +
    geom_segment(aes(x = from_V1, xend = to_V1, y = from_V2, yend = to_V2), prediction$dimred_trajectory_segments %>% as.data.frame) +
    scale_colour_manual(values = prediction$col_ann) +
    theme(legend.position = "none")
  process_dynplot(g, prediction$id)
}


