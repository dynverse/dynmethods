#' Description for StemID
#' @export
description_stemid <- function() create_description(
  name = "StemID",
  short_name = "StemID",
  package_loaded = c(),
  package_required = c("StemID"),
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
    makeIntegerParam(id = "outminc", lower = 0L, default = 5L, upper = 100L),
    makeIntegerParam(id = "outlg", lower = 0L, default = 2L, upper = 100L),
    makeNumericParam(id = "probthr", lower = -10, default = -3, upper = -1, trafo = function(x) 10^x),
    makeNumericParam(id = "thr_lower", lower = -100, default = -40, upper = -1),
    makeNumericParam(id = "thr_upper", lower = -100, default = -40, upper = -1),
    makeNumericParam(id = "outdistquant", lower = 0, default = .95, upper = 1),
    makeLogicalParam(id = "nmode", default = FALSE),
    makeNumericParam(id = "pdishuf", lower = log(100), default = log(500), upper = log(10000), trafo = function(x) round(exp(x))), # orig 2000
    makeNumericParam(id = "pthr", lower = -4, default = -2, upper = 0, trafo = function(x) 10^x),
    makeNumericParam(id = "pethr", lower = -4, default = -2, upper = 0, trafo = function(x) 10^x),
    forbidden = quote(thr_lower > thr_upper)
  ),
  properties = c(),
  run_fun = run_stemid,
  plot_fun = plot_stemid
)

run_stemid <- function(
  expression,
  clustnr = 30,
  bootnr = 50,
  metric = "pearson",
  num_cluster_method = "sat",
  SE.method = "Tibs2001SEmax",
  SE.factor = .25,
  B.gap = 50,
  cln = 0,
  FUNcluster = "kmedoids",
  dimred_method = "tsne", # tsne, sammon, tsne_initcmd
  outminc = 5,
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
  requireNamespace("StemID")

  # initialize SCseq object with transcript expression
  sc <- StemID::SCseq(data.frame(expression, check.names = FALSE))

  # filtering of expression data
  sc <- StemID::filterdata(sc, mintotal = 1, minexpr = 0, minnumber = 0, maxexpr = Inf, downsample = TRUE, dsn = 1)

  # k-medoids clustering
  do_gap <- num_cluster_method == "gap"
  do_sat <- num_cluster_method == "sat"
  sc <- StemID::clustexp(sc, clustnr = clustnr, bootnr = bootnr, metric = metric, do.gap = do_gap,
                 sat = do_sat, SE.method = SE.method, SE.factor = SE.factor,
                 B.gap = B.gap, cln = cln, FUNcluster = FUNcluster)

  # compute t-SNE map
  sammonmap <- dimred_method == "sammon"
  initial_cmd <- dimred_method == "tsne_initcmd"
  sc <- StemID::comptsne(sc, sammonmap = sammonmap, initial_cmd = initial_cmd)

  # detect outliers and redefine clusters
  thr <- 2^(thr_lower:thr_upper)
  sc <- StemID::findoutliers(
    sc,
    outminc = outminc,
    outlg = outlg,
    probthr = probthr,
    thr = thr,
    outdistquant = outdistquant
  )

  # initialization
  ltr <- StemID::Ltree(sc)

  # computation of the entropy
  ltr <- StemID::compentropy(ltr)

  # computation of the projections for all cells
  ltr <- StemID::projcells(ltr, cthr = 0, nmode = nmode)

  # computation of the projections for all cells after randomization
  ltr <- StemID::projback(ltr, pdishuf = pdishuf, nmode = nmode)

  # assembly of the lineage tree
  ltr <- StemID::lineagetree(ltr, pthr = pthr, nmode = nmode)

  # compute a spanning tree
  ltr <- StemID::compspantree(ltr)

  # get network info
  cluster_network <- data_frame(
    from = as.character(ltr@ldata$m[-1]),
    to = as.character(ltr@trl$trl$kid),
    length = ltr@trl$dc[cbind(from, to)],
    directed = FALSE
  )

  # project cells onto segments
  out <- project_cells_to_segments(
    cluster_network = cluster_network,
    cluster_space = ltr@ldata$cnl,
    sample_space = ltr@ltcoord,
    sample_cluster = as.character(ltr@ldata$lp),
    num_segments_per_edge = 100,
    milestone_rename_fun = function(x) paste0("M", x)
  )

  # get colours
  col_ann <- setNames(ltr@sc@fcol, out$milestone_ids)

  # return output
  wrap_ti_prediction(
    trajectory_type = "multifurcating",
    id = "StemID",
    cell_ids = rownames(expression),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    space = out$space_df,
    centers = out$centers_df,
    edge = out$edge_df,
    col_ann = col_ann
  )
}

plot_stemid <- function(prediction) {
  g <- ggplot() +
    geom_point(aes(V1, V2), prediction$space, size = 2, colour = "darkgray", na.rm = TRUE) +
    geom_text(aes(V1, V2, label = label, colour = label), prediction$space, na.rm = TRUE) +
    geom_text(aes(V1, V2, label = clus_id), prediction$centers, size = 8) +
    geom_segment(aes(x = from.V1, xend = to.V1, y = from.V2, yend = to.V2), prediction$edge) +
    scale_colour_manual(values = prediction$col_ann) +
    theme(legend.position = "none")
  process_dynplot(g, prediction$id)
}


