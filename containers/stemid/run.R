library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

source("StemID/RaceID2_StemID_class.R")
#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_genes = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/stemid/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   purrr::map(~ .$default)
#' params$pdishuf <- 50

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
  downsample = TRUE,
  dsn = 1
)

# k-medoids clustering
do_gap <- num_cluster_method == "gap"
do_sat <- num_cluster_method == "sat"
sc <- sc %>% clustexp(
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
sc <- sc %>% comptsne(
  sammonmap = sammonmap,
  initial_cmd = initial_cmd
)

# detect outliers and redefine clusters
sc <- sc %>% findoutliers(
  outminc = outminc,
  outlg = outlg,
  probthr = probthr,
  thr = 10^(thr_lower:thr_upper),
  outdistquant = outdistquant
)

# initialization
ltr <- Ltree(sc)

# computation of the entropy
ltr <- ltr %>% compentropy()

# computation of the projections for all cells
ltr <- ltr %>% projcells(cthr = 0, nmode = nmode)

# computation of the projections for all cells after randomization
ltr <- ltr %>% projback(pdishuf = pdishuf, nmode = nmode)

# assembly of the lineage tree
ltr <- ltr %>% lineagetree(pthr = pthr, nmode = nmode)

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
    length =  dc[cbind(from, to)],
    directed = FALSE
  ) %>%
  dplyr::select(from, to, length, directed)

output <- lst(
  milestone_ids = rownames(ltr@ldata$cnl %>% as.matrix),
  milestone_network = cluster_network,
  dimred_milestones = ltr@ldata$cnl %>% as.matrix,
  dimred = ltr@ltcoord,
  milestone_assignment_cells = as.character(ltr@ldata$lp) %>% setNames(rownames(counts)),
  timings = checkpoints
)


#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
