library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(Mpath)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (counts, groups_id, distMethod = "euclidean", method = "kmeans", 
    numcluster = 11, numcluster_null = TRUE, diversity_cut = 0.6, 
    size_cut = 0.05) 
{
    requireNamespace("igraph")
    if (numcluster_null) {
        numcluster <- NULL
    }
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    sample_info <- groups_id %>% rename(GroupID = group_id) %>% 
        as.data.frame
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    landmark_cluster <- Mpath::landmark_designation(rpkmFile = t(counts), 
        baseName = NULL, sampleFile = sample_info, distMethod = distMethod, 
        method = method, numcluster = numcluster, diversity_cut = diversity_cut, 
        size_cut = size_cut, saveRes = FALSE) %>% mutate_if(is.factor, 
        as.character)
    milestone_ids <- unique(landmark_cluster$landmark_cluster)
    if (length(milestone_ids) == 1) {
        stop("Mpath only detected one landmark")
    }
    network <- Mpath::build_network(exprs = t(counts), baseName = NULL, 
        landmark_cluster = landmark_cluster, distMethod = distMethod, 
        writeRes = FALSE)
    trimmed_network <- Mpath::trim_net(nb12 = network, writeRes = FALSE)
    class(trimmed_network) <- NULL
    milestone_network <- trimmed_network %>% reshape2::melt(varnames = c("from", 
        "to"), value.name = "length") %>% mutate_if(is.factor, 
        as.character) %>% filter(length > 0, from < to) %>% mutate(directed = FALSE)
    grouping <- with(landmark_cluster, setNames(landmark_cluster, 
        cell))
    wrap_prediction_model(cell_ids = rownames(counts), groups_id = groups_id) %>% 
        add_grouping(grouping = grouping) %>% add_cluster_graph(milestone_network = milestone_network) %>% 
        add_timings(timings = tl %>% add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')