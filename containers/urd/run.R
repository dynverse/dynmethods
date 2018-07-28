library(jsonlite)
library(readr)
library(dplyr)
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

# just laod the data, filtering has already been done by dynnormalizer
urd <- createURD(count.data = t(count), min.cells = 0, min.counts = 0, min.genes = 0)

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# Calculate Diffusion Map
urd <- calcDM(urd, knn = 100, sigma = 16)

# Then we run 'flood' simulations
urd.floods <- floodPseudotime(urd, root.cells = start_id, n = 20, minimum.cells.flooded = 2, verbose = TRUE)

# The we process the simulations into a pseudotime
urd <- floodPseudotimeProcess(urd, urd.floods, floods.name = "pseudotime")

urd@pseudotime
# all NA's?

# Calculate PCA and tSNE
urd <- calcPCA(urd, mp.factor = 1.5)
urd <- calcTsne(urd)

# Calculate graph clustering of these cells
num_nns <- 20
urd <- graphClustering(urd, num.nn = num_nns, do.jaccard=T, method="Louvain")
cluster_name <- paste0("Louvain-", num_nns)

# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
urd.ptlogistic <- pseudotimeDetermineLogistic(urd, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)

# Bias the transition matrix acording to pseudotime
urd.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(urd, "pseudotime", logistic.params=urd.ptlogistic))

# Simulate the biased random walks from each tip
urd.walks <- simulateRandomWalksFromTips(urd, tip.group.id=cluster_name, root.cells=root.cells, transition.matrix = urd.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)
# Load the cells used for each tip into the URD object
urd.tree <- loadTipCells(urd, "tip.clusters")

# Build the tree
urd.tree <- buildTree(urd.tree, pseudotime = "pseudotime", tips.use=1:2, divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)

# Name the segments based on our previous determination of the identity of tips 1 and 2.
urd.tree <- nameSegments(urd.tree, segments=c("1","2"), segment.names = c("Notochord", "Prechordal Plate"), short.names = c("Noto", "PCP"))

plotTree(urd.tree, "stage", title="Developmental Stage")
plotTree(urd.tree, "GSC", title="GSC (early prechordal plate marker)")
plotTree(urd.tree, "NOTO", title="NOTO (early notochord marker)")
plotTree(urd.tree, "HE1A", title="HE1A (prechordal plate differentiation marker")
plotTree(urd.tree, "COL8A1A", title="COL8A1A (notochord differentiation marker")


# Process the biased random walks into visitation frequencies
urd <- processRandomWalksFromTips(urd, urd.walks, verbose = F)




# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())



# return output
model <- lst(
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
