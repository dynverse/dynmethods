#' Inferring trajectories with StemID
#'
#' @inherit ti_identity description
#' @inheritParams ti_stemid2
#'
#' @export
ti_stemid <- create_ti_method(
  name = "StemID",
  short_name = "stemid",
  package_loaded = c(),
  package_required = c("StemID"),
  doi = "10.1016/j.stem.2016.05.010",
  trajectory_types = c("linear", "bifurcation", "convergence", "multifurcation"),
  topology_inference = "free",
  type = "algorithm",
  authors = list(
    list(
      given = "Dominic",
      family = "Gr\u00fcn",
      email = "gruen@ie-freibug.mpg.de",
      github = "dgrun"
    ),
    list(
      given = "Alexander",
      family = "van Oudenaarden",
      email = "a.vanoudenaarden@hubrecht.eu",
      github = "avolab"
    )
  ),
  publication_date = "2016-06-21",
  code_url = "https://github.com/dgrun/StemID",
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

    pvalue_cutoff = list(
      type = "numeric",
      default = 0.05,
      lower = 0,
      upper = 1,
      description = "Upper cutoff of p-value for determining significant edges"
    ),

    linkscore_cutoff = list(
      type = "numeric",
      default = 0.2,
      lower = 0,
      upper = 1,
      description = "Lower cutoff of the linkscore for determining significant edges"
    ),

    forbidden = "thr_lower > thr_upper"
  ),
  run_fun = "dynmethods::run_stemid",
  plot_fun = "dynmethods::plot_stemid"
)

run_stemid <- function(
  counts,
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
  pethr = .01,
  pvalue_cutoff,
  linkscore_cutoff
) {
  requireNamespace("StemID")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # initialize SCseq object with transcript expression
  sc <- StemID::SCseq(data.frame(t(counts), check.names = FALSE))

  # filtering of expression data
  sc <- sc %>% StemID::filterdata(
    mintotal = 1,
    minexpr = 0,
    minnumber = 0,
    maxexpr = Inf,
    downsample = TRUE,
    dsn = 1
  )

  # k-medoids clustering
  do_gap <- num_cluster_method == "gap"
  do_sat <- num_cluster_method == "sat"
  sc <- sc %>% StemID::clustexp(
    clustnr = clustnr,
    bootnr = bootnr,
    metric = metric,
    do.gap = do_gap,
    sat = do_sat,
    SE.method = SE.method,
    SE.factor = SE.factor,
    B.gap = B.gap,
    cln = cln,
    FUNcluster = FUNcluster
  )

  # compute t-SNE map
  sammonmap <- dimred_method == "sammon"
  initial_cmd <- dimred_method == "tsne_initcmd"
  sc <- sc %>% StemID::comptsne(
    sammonmap = sammonmap,
    initial_cmd = initial_cmd
  )

  # detect outliers and redefine clusters
  sc <- sc %>% StemID::findoutliers(
    outminc = outminc,
    outlg = outlg,
    probthr = probthr,
    thr = 10^(thr_lower:thr_upper),
    outdistquant = outdistquant
  )

  # initialization
  ltr <- StemID::Ltree(sc)

  # computation of the entropy
  ltr <- ltr %>% StemID::compentropy()

  # computation of the projections for all cells
  ltr <- ltr %>% StemID::projcells(cthr = 0, nmode = nmode)

  # computation of the projections for all cells after randomization
  ltr <- ltr %>% StemID::projback(pdishuf = pdishuf, nmode = nmode)

  # assembly of the lineage tree
  ltr <- ltr %>% StemID::lineagetree(pthr = pthr, nmode = nmode)

  # compute a spanning tree
  ltr <- ltr %>% StemID::compspantree()

  # compute p value
  ltr <- ltr %>% StemID:::comppvalue()

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # get linkscores and pvalues
  cluster_network_linkscore <- ltr@cdata$linkscore %>%
    tibble::rownames_to_column("from") %>%
    tidyr::gather("to", "linkscore", -from)

  cluster_network_pvalue <- ltr@cdata$pvn.e %>%
    tibble::rownames_to_column("from") %>%
    tidyr::gather("to", "pvalue", -from)

  # combine into one cluster network
  cluster_network <- left_join(
    cluster_network_linkscore,
    cluster_network_pvalue,
    c("from", "to")
  ) %>%
    mutate_at(c("from", "to"), ~gsub("cl\\.(.*)", "\\1", .))

  # filter the cluster network
  cluster_network <- cluster_network %>%
    filter(
      pvalue <= pvalue_cutoff,
      linkscore >= linkscore_cutoff
    )

  # get distances between clusters
  cluster_network <- cluster_network %>%
    mutate(
      length =  ltr@trl$dc[cbind(from, to)],
      directed = FALSE
    ) %>%
    select(from, to, length, directed)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_dimred_projection(
    milestone_ids = rownames(ltr@ldata$cnl %>% as.matrix),
    milestone_network = cluster_network,
    dimred_milestones = ltr@ldata$cnl %>% as.matrix,
    dimred = ltr@ltcoord,
    milestone_assignment_cells = as.character(ltr@ldata$lp) %>% setNames(rownames(counts)),
    num_segments_per_edge = 100,
    col_ann = ltr@sc@fcol
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_stemid <- function(prediction) {
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


