library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(TSCAN)
library(igraph)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (counts, minexpr_percent = 0, minexpr_value = 0, cvcutoff = 0, 
    clusternum_lower = 2, clusternum_upper = 9, modelNames = "VVV") 
{
    requireNamespace("TSCAN")
    requireNamespace("igraph")
    clusternum <- seq(clusternum_lower, clusternum_upper, 1)
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    cds_prep <- TSCAN::preprocess(t(as.matrix(counts)), takelog = TRUE, 
        logbase = 2, pseudocount = 1, clusternum = NULL, minexpr_value = minexpr_value, 
        minexpr_percent = minexpr_percent, cvcutoff = cvcutoff)
    cds_clus <- TSCAN::exprmclust(cds_prep, clusternum = clusternum, 
        modelNames = modelNames, reduce = TRUE)
    cds_order <- TSCAN::TSCANorder(cds_clus)
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    cluster_network <- cds_clus$MSTtree %>% igraph::as_data_frame() %>% 
        rename(length = weight) %>% mutate(directed = FALSE)
    sample_space <- cds_clus$pcareduceres
    cluster_space <- cds_clus$clucenter
    rownames(cluster_space) <- as.character(seq_len(nrow(cluster_space)))
    colnames(cluster_space) <- colnames(sample_space)
    wrap_prediction_model(cell_ids = rownames(counts)) %>% add_dimred_projection(milestone_ids = rownames(cluster_space), 
        milestone_network = cluster_network, dimred_milestones = cluster_space, 
        dimred = sample_space, milestone_assignment_cells = cds_clus$clusterid, 
        num_segments_per_edge = 100) %>% add_timings(timings = tl %>% 
        add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')