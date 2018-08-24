#!/usr/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(gng)
library(igraph)
library(dyndimred)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

expression <- data$expression
#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# perform dimensionality reduction
space <- dyndimred::dimred(expression, method = params$dimred, ndim = params$ndim)

# calculate GNG
gng_out <- gng::gng(
  space,
  max_iter = params$max_iter,
  max_nodes = params$max_nodes,
  assign_cluster = FALSE
)
node_dist <- stats::dist(gng_out$node_space) %>% as.matrix

# transform to milestone network
node_names <- gng_out$nodes %>% mutate(name = as.character(name))
milestone_network <- gng_out$edges %>%
  select(from = i, to = j) %>%
  mutate(
    length = node_dist[cbind(from, to)],
    directed = FALSE
  ) %>%
  select(from, to, length, directed)

# apply MST, if so desired
if (params$apply_mst) {
  gr <- igraph::graph_from_data_frame(milestone_network, directed = F, vertices = node_names$name)
  milestone_network <- igraph::minimum.spanning.tree(gr, weights = igraph::E(gr)$length) %>% igraph::as_data_frame()
}

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# return output
output <- lst(
  cell_ids = rownames(expression),
  milestone_ids = rownames(gng_out$node_space),
  milestone_network,
  dimred_milestones = gng_out$node_space,
  dimred = space,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
