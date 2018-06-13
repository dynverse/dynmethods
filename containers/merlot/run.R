library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(merlot)
library(destiny)

checkpoints <- list()

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')
expression <- data$expression
end_n <- data$end_n
start_id <- data$start_id


#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

checkpoints$method_afterpreproc <- as.numeric(Sys.time())

# create n_local vector
n_local <- seq(params$n_local_lower, params$n_local_upper, by = 1)

#### Example fromrom inst/examples/ExampleGuo2010.R
if(!is.null(end_n)) {
  n_components_to_use <- end_n - 1
}
ndim <- max(n_components_to_use, params$ndim) # always make sure that enough components are extracted, even if the provided n_components is too low

# Embed Cells into their manifold, in this case we use Diffusion Maps as calculated by Destiny
DatasetDM <- destiny::DiffusionMap(
  data = expression,
  sigma = params$sigma,
  distance = params$distance,
  n_eigs = ndim,
  density_norm = params$density_norm,
  n_local = n_local,
  verbose = F
)

# Extract dimensionality reduction
CellCoordinates <- DatasetDM@eigenvectors[,seq_len(n_components_to_use)]

# We calculate the scaffold tree using the first 3 diffusion components from the diffusion map
ScaffoldTree <- merlot::CalculateScaffoldTree(
  CellCoordinates = CellCoordinates,
  NEndpoints = end_n
)

# Set the number of nodes to be used to build the Principal Elastic Tree.
# This is now a parameter of the method

# We calculate the elastic principal tree using the scaffold tree for its initialization
ElasticTree <- merlot::CalculateElasticTree(
  ScaffoldTree = ScaffoldTree,
  N_yk = params$N_yk,
  lambda_0 = params$lambda_0,
  mu_0 = params$mu_0,
  FixEndpoints = params$FixEndpoints
)

# Embedd the principal elastic tree into the gene expression space from which it was calculated.
EmbeddedTree <- merlot::GenesSpaceEmbedding(
  ExpressionMatrix = expression,
  ElasticTree = ElasticTree,
  lambda_0 = params$lambda_0,
  mu_0 = params$mu_0,
  increaseFactor_mu = params$increaseFactor_mu,
  increaseFactor_lambda = params$increaseFactor_lambda
)

# Calculate Pseudotimes for the nodes in the Tree in the full gene expression space.
# T0=3 means that the Endpoint number 3 in the Endpoints list corresponds to the zygote fate and is used as initial pseudotime t0
# Any given cell can be used as t0 by specifying its index using the parameter C0=cell_index
if (is.null(start_id)) {
  Pseudotimes <- merlot::CalculatePseudotimes(EmbeddedTree, T0=1)
} else {
  Pseudotimes <- merlot::CalculatePseudotimes(EmbeddedTree, C0=which(rownames(expression) == start_id))
}

checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Save output & process output                                            ####

# first add both the milestone network (without lengths) and progressions
milestone_network <- ElasticTree$Edges %>%
  as.data.frame() %>%
  purrr::set_names(c("from", "to")) %>%
  mutate_all(as.character) %>%
  mutate(edge_id = row_number())

progressions <- tibble(
  cell_id = rownames(expression),
  edge_id = ElasticTree$Cells2Branches,
  pseudotime = Pseudotimes$Proyected_Times_Cells
) %>%
  left_join(milestone_network, "edge_id")

# now calculate milestone network lengths
milestone_network <- left_join(
  milestone_network,
  progressions %>%
    group_by(edge_id) %>%
    summarise(length = max(pseudotime) - min(pseudotime)),
  "edge_id"
) %>%
  mutate(length = ifelse(is.na(length), mean(length, na.rm=T), length)) %>%
  mutate(directed = TRUE) %>%
  select(from, to, length, directed)

# now calculate percentages of progression
progressions <- progressions %>%
  group_by(edge_id) %>%
  mutate(percentage = (pseudotime - min(pseudotime))/(max(pseudotime) - min(pseudotime))) %>%
  ungroup() %>%
  select(cell_id, from, to, percentage)

# get dimred
dimred <- CellCoordinates %>%
  as.data.frame()
dimred$cell_id <- rownames(expression)

# save
write_csv(milestone_network, '/output/milestone_network.csv')
write_csv(progressions, '/output/progressions.csv')
write_csv(dimred, '/output/dimred.csv')
