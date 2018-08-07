library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

source("/RaceID3_StemID2/RaceID3_StemID2_class.R")

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/stemid2/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts
list2env(params, .GlobalEnv)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# initialize SCseq object with transcript expression
sc <- SCseq(data.frame(t(counts), check.names = FALSE))

# filtering of expression data
sc <- sc %>% filterdata(
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
  sc %>% clustexp(
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
  sc %>% clustexp(
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
sc <- sc %>% comptsne(
  sammonmap = sammonmap,
  initial_cmd = initial_cmd,
  fast = TRUE, # newly added
  perplexity = 30 # newly added
)

# detect outliers and redefine clusters
sc <- sc %>% findoutliers(
  outminc = 5,
  outlg = outlg,
  probthr = probthr,
  thr = 10^(thr_lower:thr_upper),
  outdistquant = outdistquant
)

sc <- sc %>% rfcorrect(
  final = TRUE, # newly added
  nbfactor = 5 # newly added
)

# initialization
ltr <- Ltree(sc)

# computation of the entropy
ltr <- ltr %>% compentropy()

# computation of the projections for all cells
ltr <- ltr %>% projcells(
  cthr = 0, # default = 2
  nmode = nmode
)

# computation of the projections for all cells after randomization
ltr <- ltr %>% projback(
  pdishuf = pdishuf,
  nmode = nmode,
  fast = FALSE # newly added
)

# assembly of the lineage tree
ltr <- ltr %>% lineagetree(
  pthr = pthr,
  nmode = nmode,
  fast = FALSE # newly added
)

# compute a spanning tree
medoids <- compmedoids(ltr@sc@fdata, ltr@sc@cpart)
cent <- ltr@sc@fdata[,medoids]
dc <- as.data.frame(1 - cor(cent))
names(dc) <- sort(unique(ltr@sc@cpart))
rownames(dc) <- sort(unique(ltr@sc@cpart))

# compute p value
ltr <- ltr %>% comppvalue()

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())


# get linkscores and pvalues
milestone_network_linkscore <- ltr@cdata$linkscore %>%
  tibble::rownames_to_column("from") %>%
  tidyr::gather("to", "linkscore", -from)

milestone_network_pvalue <- ltr@cdata$pvn.e %>%
  tibble::rownames_to_column("from") %>%
  tidyr::gather("to", "pvalue", -from)

# combine into one cluster network
milestone_network <- left_join(
  milestone_network_linkscore,
  milestone_network_pvalue,
  c("from", "to")
) %>%
  mutate_at(c("from", "to"), ~gsub("cl\\.(.*)", "\\1", .)) %>%
  filter(
    pvalue <= pvalue_cutoff,
    linkscore >= linkscore_cutoff
  ) %>%
  mutate(
    length =  dc[cbind(from, to)],
    directed = FALSE
  ) %>%
  dplyr::select(from, to, length, directed)


dimred_milestones <- ltr@ldata$cnl %>% as.matrix
dimred <- ltr@ltcoord %>% na.omit
cell_ids <- rownames(dimred)
milestone_ids <- rownames(dimred_milestones)
grouping <- ltr@ldata$lp[cell_ids]

output <- lst(
  cell_ids,
  milestone_ids,
  milestone_network = milestone_network,
  dimred_milestones,
  dimred,
  grouping,
  timings = checkpoints
)

#' @examples
#' dynwrap::wrap_data(
#'   cell_ids = cell_ids
#' ) %>% dynwrap::add_dimred_projection(
#'   milestone_ids = milestone_ids,
#'   milestone_network = milestone_network,
#'   dimred_milestones = dimred_milestones,
#'   dimred = dimred,
#'   grouping = grouping
#' )

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
