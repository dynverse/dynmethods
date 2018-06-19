#' Inferring trajectories with StemID2
#'
#' @param clustnr maximum number of clusters for the computation of the gap statistic or derivation of the cluster number by the saturation criterion. Default is 30. If more major cell types are expected a higher number should be chosen.
#' @param bootnr number of booststrapping runs for clusterboot. Default is 50
#' @param metric the input data are transformed to a distance object. Distances can be computed based on different metrics. Possible values are "pearson", "spearman", "logpearson", "euclidean", "kendall", "maximum", "manhattan", "canberra", "binary" or "minkowski". Default is "pearson". In case of the correlation based methods, the distance is computed as 1 - correlation. K-medoids clustering is performed on this distance object.
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
ti_stemid2 <- create_ti_method(
  name = "StemID2",
  short_name = "stemid2",
  package_loaded = c(),
  package_required = c("StemID2"),
  doi = "10.1016/j.stem.2016.05.010",
  trajectory_types = c("linear", "bifurcation", "convergence", "multifurcation", "binary_tree", "tree"),
  topology_inference = "free",
  type = "algorithm",
  authors = list(
    list(
      given = "Dominic",
      family = "Grün",
      email = "gruen@ie-freibug.mpg.de",
      github = "dgrun"
    ),
    list(
      given = "Alexander",
      family = "van Oudenaarden",
      email = "a.vanoudenaarden@hubrecht.eu"
    )
  ),
  code_url = "https://github.com/dgrun/RaceID3_StemID2",
  parameters = list(
    clustnr = list(
      type = "integer",
      default = 30L,
      upper = 100L,
      lower = 10L,
      description = "maximum number of clusters for the computation of the gap statistic or derivation of the cluster number by the saturation criterion. Default is 30. If more major cell types are expected a higher number should be chosen."),

    bootnr = list(
      type = "integer",
      default = 50L,
      upper = 100L,
      lower = 20L,
      description = "number of booststrapping runs for clusterboot. Default is 50"),

    metric = list(
      type = "discrete",
      default = "pearson",
      values = c("pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
      description = "the input data are transformed to a distance object. Distances can be computed based on different metrics. Possible values are \"pearson\", \"spearman\", \"logpearson\", \"euclidean\", \"kendall\", \"maximum\", \"manhattan\", \"canberra\", \"binary\" or \"minkowski\". Default is \"pearson\". In case of the correlation based methods, the distance is computed as 1 - correlation. K-medoids clustering is performed on this distance object."),

    num_cluster_method = list(
      type = "discrete",
      default = "sat",
      values = c("sat", "gap", "manual"),
      description = "the type of clustering method, can be sat, gap or manual"),
    SE.method = list(
      type = "discrete",
      default = "Tibs2001SEmax",
      values = c("firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax"),
      description = "the clustering routine calls a modified version of the maxSE function from the cluster package to determine the first local maximum of the gap statistic. By default, we use the method \"Tibs2001SEmax\" for calling the first local maximum (see specification of maxSE). This method requires that the maximum exceeds the values of its neighbors by a fraction of their standard deviation. This fraction is defined by the parameter SE.factor. All methods defined for the original maxSE function can also be used."),

    SE.factor = list(
      type = "numeric",
      default = 0.25,
      upper = 1,
      lower = 0,
      description = "fraction of the standard deviation by which the local maximum is required to differ from the neighboring points it is compared to. Default is 0.25."),
    B.gap = list(
      type = "integer",
      default = 50L,
      upper = 100L,
      lower = 20L,
      description = "number of bootstrap runs for the calculation of the gap statistic. Default is 50."),
    cln = list(
      type = "integer",
      default = 30L,
      upper = 100L,
      lower = 20L,
      description = "the number of clusters for k-medoids clustering. Default is 0. In this case, the cluster number is determined based on the gap statistic6 and do.gap has to be TRUE."),

    FUNcluster = list(
      type = "discrete",
      default = "kmedoids",
      values = c("kmedoids", "kmeans", "hclust"),
      description = "the clustering method applied. One of the following methods can be selected: kmedoids, kmeans, hclust. RaceID3 is designed for k-medoids clustering and therefore it is recommended to use only the kmedoids method. Default is kmedoids."),
    dimred_method = list(
      type = "discrete",
      default = "tsne",
      values = c("tsne", "sammon", "tsne_initcmd"),
      description = "the dimensionality reduction method, can be tsne, sammon or tsne_initcmd"),

    outminc = list(
      type = "integer",
      default = 0L,
      upper = 100L,
      lower = 0L,
      description = "expression cutoff for the identification of outlier genes is defined. Default is 5."),
    outlg = list(
      type = "integer",
      default = 2L,
      upper = 100L,
      lower = 0L,
      description = "minimal number of outlier genes required to identify a cell as an outlier. Default is 2."),
    probthr = list(
      type = "numeric",
      default = 1e-3,
      upper = 1e-1,
      lower = 1e-10,
      description = "defines the probability threshold for outlier calling. If the probability of observing a given expression level for a gene in a cell is lower than this cutoff (based on the negative binomial distribution for the calibrated noise model), the cell is considered an outlier for this gene. Default is 10-3."),
    thr_lower = list(
      type = "integer",
      default = -10,
      upper = -1,
      lower = -100,
      description = "lower probability for which the number of outliers is computed in order to plot the dependence of the number of outliers on the probability threshold"),
    thr_upper = list(
      type = "numeric",
      default = -5,
      upper = -1,
      lower = -100,
      description = "upper probability for which the number of outliers is computed in order to plot the dependence of the number of outliers on the probability threshold"),
    outdistquant = list(

      type = "numeric",
      default = 0.95,
      upper = 1,
      lower = 0,
      description = "outlier cells are merged to outlier clusters if their similarity exceeds the outdistquant-quantile of the similarity distribution for all pairs of cells that are together in one of the original clusters. Default is 0.95."),

    nmode = list(
      type = "logical",
      default = FALSE,
      values = c("TRUE", "FALSE"),
      description = "Boolean argument. If nmode is set to TRUE the assignment to inter-cluster links for each cell is not done based on the longest projection, but based on identifying the cluster (other than the cluster the cell belongs to) that contains the nearest neighbor of the cell, i. e. the cell with the most similar transcriptome. The coordinate on the assigned link is still derived based on the projection. Default is FALSE."),

    pdishuf = list(
      type = "integer",
      default = 2000,
      upper = 10000,
      lower = 100,
      description = "positive integer. This is the number of randomizations to be performed. As a rule of thumb this number should be at least one order of magnitude larger than the desired p-value on the significance of the number of cells on a connection. Default is 2000."),
    pthr = list(
      type = "numeric",
      default = 1e-2,
      upper = 1e-4,
      lower = 0,
      description = "positive number. This number corresponds to the p-value threshold, which is used to determine, whether the magnitude of an observed trajectory is significantly larger than observed for the randomized background distribution. This criterion is not used to infer significance of a link, but shown in a graphical representation of the tree"),

    pethr = list(
      type = "numeric",
      default = 1e-2,
      upper = 1e-4,
      lower = 0,
      description = "positive number. This number corresponds to the p-value threshold, which isused to determine for each link if it is populated by a number of cellssignificantly larger than expected for the randomized background distribution. This p-value threshold determines, which connections are considered validdifferentiation trajectories in the derived lineage tree."),
    forbidden = "thr_lower > thr_upper"
  ),
  run_fun = "dynmethods::run_stemid2",
  plot_fun = "dynmethods::plot_stemid2"
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
    thr = 10^(thr_lower:thr_upper),
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


