library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(StemID)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  counts,
  clustnr = 30,
  bootnr = 50,
  metric = "pearson",
  num_cluster_method = "sat",
  SE.method = "Tibs2001SEmax",
  SE.factor = 0.25,
  B.gap = 50,
  cln = 30,
  FUNcluster = "kmedoids",
  dimred_method = "tsne",
  outminc = 0,
  outlg = 2,
  probthr = 0.001,
  thr_lower = -10,
  thr_upper = -5,
  outdistquant = 0.95,
  nmode = FALSE,
  pdishuf = 2000,
  pthr = 1e-04,
  pethr = 1e-04,
  pvalue_cutoff = 0.05,
  linkscore_cutoff = 0.2
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

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')