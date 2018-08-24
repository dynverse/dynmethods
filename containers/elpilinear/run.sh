#!/usr/bin/Rscript

library(ElPiGraph.R)
library(dplyr)
library(purrr)
library(readr)


#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "binary_tree") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/elpilinear/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

#   ____________________________________________________________________________
#   Infer the trajectory                                                    ####
# choose graph topology
if (!is.null(data$trajectory_type)) {
  principal_graph_function <- switch(
    data$trajectory_type,
    linear = computeElasticPrincipalCurve,
    directed_linear = computeElasticPrincipalCurve,
    undirected_linear = computeElasticPrincipalCurve,
    cycle = computeElasticPrincipalCircle,
    directed_cycle = computeElasticPrincipalCircle,
    undirected_cycle = computeElasticPrincipalCircle,
    computeElasticPrincipalTree
  )
} else {
  principal_graph_function <- switch(
    params$topology,
    linear = computeElasticPrincipalCurve,
    cycle = computeElasticPrincipalCircle,
    computeElasticPrincipalTree
  )
}

# infer the principal graph, from https://github.com/Albluca/ElPiGraph.R/blob/master/guides/base.md
principal_graph <- principal_graph_function(
  X = expression,
  NumNodes = params$NumNodes,
  NumEdges = params$NumEdges,
  InitNodes = params$InitNodes,
  MaxNumberOfIterations = params$MaxNumberOfIterations,
  eps = params$eps,
  CenterData = params$CenterData,
  Lambda = params$Lambda,
  Mu = params$Mu,
  drawAccuracyComplexity = FALSE,
  drawEnergy = FALSE,
  drawPCAView = FALSE,
  n.cores = 1
)

# compute pseudotime, from https://github.com/Albluca/ElPiGraph.R/blob/master/guides/pseudo.md
PartStruct <- PartitionData(
  X = expression,
  NodePositions = principal_graph[[1]]$NodePositions,
  nCores = 1
)

ProjStruct <- project_point_onto_graph(
  X = expression,
  NodePositions = principal_graph[[1]]$NodePositions,
  Edges = principal_graph[[1]]$Edges$Edges,
  Partition = PartStruct$Partition
)

checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Process & save the model                                               ####
milestone_network <- ProjStruct$Edges %>%
  as_data_frame() %>%
  rename(from = row, to = col) %>%
  mutate(
    from = paste0("M", from),
    to = paste0("M", to),
    length = ProjStruct$EdgeLen,
    directed = FALSE
  )

progressions <- tibble(cell_id = rownames(expression), edge_id = ProjStruct$EdgeID) %>%
  left_join(milestone_network %>% select(from, to) %>% mutate(edge_id = row_number()), "edge_id") %>%
  select(-edge_id) %>%
  mutate(percentage = pmin(1, pmax(0, ProjStruct$ProjectionValues)))

output <- lst(
  cell_ids = rownames(expression),
  milestone_network,
  progressions,
  timings = checkpoints
)

write_rds(output, "/ti/output/output.rds")
