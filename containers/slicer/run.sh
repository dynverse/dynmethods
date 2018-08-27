#!/usr/local/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(SLICER)
library(lle)
library(igraph)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")


#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/slicer/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression
start_id <- data$start_id
features_id <- data$features_id
end_id <- data$end_id

start_cell <- sample(start_id, 1)

if (!is.null(features_id)) {
  # use 'neighbourhood variance' to identify genes that vary smoothly
  features_id <- SLICER::select_genes(expression)
}
expr_filt <- expression[, features_id]

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# determine k for knn
k <- SLICER::select_k(expr_filt, kmin = params$kmin)

# perform local linear embedding
traj_lle <- lle::lle(expr_filt, m = params$m, k = params$k)$Y
rownames(traj_lle) <- rownames(expr_filt)
colnames(traj_lle) <- paste0("comp_", seq_len(ncol(traj_lle)))

# get LLE KNN graph
traj_graph <- SLICER::conn_knn_graph(traj_lle, k = params$k)
igraph::V(traj_graph)$name <- rownames(traj_lle)

# find extreme cells
if (is.null(end_id)) {
  ends <- SLICER::find_extreme_cells(traj_graph, traj_lle)
} else {
  ends <- match(c(start_cell, end_id), rownames(expression))
}

# order cells
start <- which(rownames(expr_filt) == start_cell)
cells_ordered <- SLICER::cell_order(traj_graph, start)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# get shortest paths to start and all other nodes
shortest_paths <- igraph::shortest_paths(traj_graph, start)
gr <- lapply(shortest_paths$vpath, function(path) {
  P <- rbind(path[-length(path)], path[-1]) %>% as.vector
  igraph::E(traj_graph, P = P)
})
subgr <- igraph::subgraph.edges(traj_graph, eids = unique(unlist(gr)))

# prepare sample graph simplification
cell_graph <- igraph::as_data_frame(subgr, "edges") %>%
  dplyr::select(from, to, length = weight) %>%
  mutate(
    directed = TRUE
  )
sh_p_to_ends <- igraph::shortest_paths(subgr, start, ends)
nodes_to_keep <- unique(sh_p_to_ends$vpath %>% unlist)
cell_ids <- names(igraph::V(traj_graph))
to_keep <- setNames(seq_along(cell_ids) %in% nodes_to_keep, cell_ids)

# return output
output <- lst(
  cell_ids,
  cell_graph,
  to_keep,
  dimred = traj_lle,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
