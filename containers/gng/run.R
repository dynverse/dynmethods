library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(GNG)
library(igraph)
library(dyndimred)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, dimred = "pca", ndim = 5, max_iter = 13.8155105579643, 
    max_nodes = 8, apply_mst = TRUE) 
{
    requireNamespace("GNG")
    requireNamespace("igraph")
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    space <- dyndimred::dimred(expression, method = dimred, ndim = ndim)
    gng_out <- GNG::gng(space, max_iter = max_iter, max_nodes = max_nodes, 
        assign_cluster = FALSE)
    node_dist <- stats::dist(gng_out$node_space) %>% as.matrix
    node_names <- gng_out$nodes %>% mutate(name = as.character(name))
    milestone_network <- gng_out$edges %>% select(from = i, to = j) %>% 
        mutate(length = node_dist[cbind(from, to)], directed = FALSE) %>% 
        select(from, to, length, directed)
    if (apply_mst) {
        gr <- igraph::graph_from_data_frame(milestone_network, 
            directed = F, vertices = node_names$name)
        milestone_network <- igraph::minimum.spanning.tree(gr, 
            weights = igraph::E(gr)$length) %>% igraph::as_data_frame()
    }
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_dimred_projection(milestone_ids = rownames(gng_out$node_space), 
            milestone_network = milestone_network, dimred_milestones = gng_out$node_space, 
            dimred = space) %>% add_timings(tl %>% add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')