library(jsonlite)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# Hotfix for drop = FALSE problem in URD
URD:::floodPseudotimeCalc %>%
  deparse() %>%
  gsub("cells.visited], 1, combine.probs)", "cells.visited, drop = FALSE], 1, combine.probs)", ., fixed = TRUE) %>%
  parse(text = .) %>%
  eval(envir = environment(URD:::floodPseudotimeCalc)) %>%
  utils::assignInNamespace("floodPseudotimeCalc", ., ns = "URD")

library(URD)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_features = 300, model = "bifurcating") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/urd/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
count <- data$count
start_id <- data$start_id

if (params$sigma.use <= 0) {
  params$sigma.use <- NULL
}
if (params$knn <= 0) {
  params$knn <- destiny:::find_dm_k(nrow(count))
  if (params$knn >= nrow(count)) {
    params$knn <- round(log10(nrow(count)) * 10)
  }
}

# just laod the data, filtering has already been done by dynnormalizer
urd <- createURD(count.data = t(count), min.cells = 0, min.counts = 0, min.genes = 0)

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# Calculate Diffusion Map
urd <- calcDM(
  urd,
  knn = params$knn,
  sigma.use = params$sigma.use,
  distance = params$distance
)

# Then we run 'flood' simulations
urd.floods <- floodPseudotime(
  urd,
  root.cells = start_id,
  n = params$n_floods,
  minimum.cells.flooded = 0,
  verbose = FALSE
)

# The we process the simulations into a pseudotime
urd <- floodPseudotimeProcess(
  urd,
  urd.floods,
  floods.name = "pseudotime",
  stability.div = params$stability.div
)

# Calculate PCA and tSNE
urd <- calcPCA(
  urd,
  mp.factor = params$mp.factor
)
urd <- calcTsne(
  urd,
  perplexity = params$perplexity,
  theta = params$theta,
  max_iter = params$max_iter
)

# Calculate graph clustering of these cells
urd <- graphClustering(
  urd,
  num.nn = params$num.nn,
  do.jaccard = params$do.jaccard,
  method = "Louvain"
)
cluster_name <- paste0("Louvain-", params$num.nn)

# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
urd.ptlogistic <- pseudotimeDetermineLogistic(
  urd,
  "pseudotime",
  optimal.cells.forward = params$optimal.cells.forward,
  max.cells.back = params$max.cells.back,
  do.plot = FALSE
)

# Bias the transition matrix acording to pseudotime
urd.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(
  urd,
  "pseudotime",
  logistic.params = urd.ptlogistic
))

# Simulate the biased random walks from each tip
urd.walks <- simulateRandomWalksFromTips(
  urd,
  tip.group.id = cluster_name,
  root.cells = start_id,
  transition.matrix = urd.biased.tm,
  n.per.tip = params$n.per.tip,
  root.visits = params$root.visits,
  max.steps = 5000,
  verbose = FALSE
)

# Process the biased random walks into visitation frequencies
urd <- processRandomWalksFromTips(
  urd,
  urd.walks,
  n.subsample = params$n.subsample,
  verbose = FALSE
)

# Load the cells used for each tip into the URD object
urd.tree <- loadTipCells(urd, cluster_name)

tips.use <- unique(urd.tree@group.ids[,cluster_name])

# Build the tree
urd.tree <- buildTree(
  urd.tree,
  pseudotime = "pseudotime",
  tips.use = tips.use,
  divergence.method = params$divergence.method,
  cells.per.pseudotime.bin = params$cells.per.pseudotime.bin,
  bins.per.pseudotime.window = params$bins.per.pseudotime.window,
  p.thresh = params$p.thresh,
  dendro.cell.jitter = 0,
  dendro.cell.dist.to.tree = 0,
  save.all.breakpoint.info = TRUE
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

tree_layout <- urd.tree@tree$tree.layout
subs <- tree_layout %>% filter(do.mean)
tree_layout <- tree_layout %>%
  filter(!do.mean) %>%
  left_join(subs %>% select(node.1 = node.2, node.new = node.1), by = "node.1") %>%
  mutate(node.1 = ifelse(is.na(node.new), node.1, node.new)) %>%
  select(from = node.1, to = node.2, x1, x2, y1, y2)

cell_layout <- urd.tree@tree$cell.layout

comb <- crossing(cell_layout, tree_layout) %>%
  filter(x1 <= x & x <= x2 & y1 <= y & y <= y2) %>%
  group_by(cell) %>%
  slice(1) %>%
  ungroup()

progressions <- comb %>%
  mutate(percentage = (y - y1) / (y2 - y1)) %>%
  select(cell_id = cell, from, to, percentage)

# collect milestone network
milestone_network <- tree_layout %>%
  mutate(
    length = abs(y1 - y2),
    directed = FALSE
  ) %>%
  select(from, to, length, directed)


#' @examples
#' library(dynwrap)
#' traj <- wrap_data(
#'   cell_ids = urd.tree@tree$cell.layout$cell
#' ) %>% add_trajectory(
#'   milestone_network = milestone_network,
#'   progressions = progressions
#' )
#' dynplot::plot_graph(traj)


# return output
output <- lst(
  milestone_network = milestone_network,
  progressions = progressions,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
