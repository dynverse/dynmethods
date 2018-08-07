library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(RaceID)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/raceid_stemid/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts

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
  knn = params$knn,
  ccor = params$ccor
)

# compute pairwise distances
sc <- sc %>% compdist(
  metric = params$metric,
  FSelect = FALSE
)

# perform clustering
sc <- sc %>% clustexp(
  sat = params$sat,
  samp = params$samp,
  cln = params$cln,
  clustnr = params$clustnr,
  bootnr = params$bootnr,
  FUNcluster = params$FUNcluster
)

# detect outliers and redefine clusters
sc <- sc %>% findoutliers(
  probthr = params$probthr,
  outminc = params$outminc,
  outlg = params$outlg,
  outdistquant = params$outdistquant
)

# compute t-SNE map
sc <- sc %>% comptsne(
  initial_cmd = params$initial_cmd,
  perplexity = params$perplexity
)

# initialization
ltr <- Ltree(sc)

# computation of the entropy
ltr <- ltr %>% compentropy()

# computation of the projections for all cells
ltr <- ltr %>% projcells(
  cthr = params$cthr,
  nmode = params$nmode,
  knn = params$projcells_knn,
  fr = params$fr
)

# computation of the projections for all cells after randomization
ltr <- ltr %>% projback(
  pdishuf = params$pdishuf,
  fast = params$fast
)

# assembly of the lineage tree
ltr <- ltr %>% lineagegraph()

# compute p-values for link significance
ltr <- ltr %>% comppvalue(
  pthr = params$pthr
)


# compute p value
ltr <- ltr %>% comppvalue()

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# collect information on dimreds
dimred_milestones <- ltr@ldata$cnl %>% as.matrix
rownames(dimred_milestones) <- paste0("M", ltr@ldata$m)
dimred <- ltr@ltcoord %>% na.omit
cell_ids <- rownames(dimred)
milestone_ids <- rownames(dimred_milestones)
grouping <- paste0("M", ltr@ldata$lp[cell_ids])

# calculate distance between milestones
dist_milestones <- as.matrix(dist(dimred_milestones))

# fetch milestone network by filtering the linkscore
milestone_network <- ltr@cdata$linkscore %>%
  as.matrix() %>%
  reshape2::melt(varnames = c("from", "to"), value.name = "linkscore") %>%
  na.omit() %>%
  filter(linkscore >= params$scthr) %>%
  mutate_at(c("from", "to"), ~gsub("cl.", "M", ., fixed = TRUE)) %>%
  mutate(
    length =  dist_milestones[cbind(from, to)],
    directed = FALSE
  ) %>%
  dplyr::select(from, to, length, directed)


output <- lst(
  cell_ids,
  milestone_ids,
  milestone_network,
  dimred_milestones,
  dimred,
  grouping,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
