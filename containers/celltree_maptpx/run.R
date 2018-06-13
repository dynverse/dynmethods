library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(cellTree)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, start_id = NULL, groups_id = NULL, method = "maptpx", 
    num_topics_lower = 2, num_topics_upper = 15, num_topics = NULL, 
    sd_filter = 0.5, tot_iter = 1e+06, tolerance = 0.05, absolute_width = 0, 
    width_scale_factor = 1.5, outlier_tolerance_factor = 0.1, 
    rooting_method = "null") 
{
    requireNamespace("cellTree")
    start_cell <- if (!is.null(start_id)) {
        sample(start_id, 1)
    }
    else {
        NULL
    }
    if (rooting_method == "null") {
        rooting_method <- NULL
    }
    if (is.null(num_topics)) {
        num_topics <- seq(num_topics_lower, num_topics_upper)
    }
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    lda_out <- cellTree::compute.lda(t(expression) + min(expression) + 
        1, k.topics = num_topics, method = method, log.scale = FALSE, 
        sd.filter = sd_filter, tot.iter = tot_iter, tol = tolerance)
    backbone_params <- list(lda.results = lda_out, absolute.width = absolute_width, 
        width.scale.factor = width_scale_factor, outlier.tolerance.factor = outlier_tolerance_factor, 
        rooting.method = rooting_method, only.mst = FALSE, merge.sequential.backbone = FALSE)
    if (!is.null(groups_id)) {
        backbone_params$grouping <- groups_id %>% slice(match(cell_id, 
            rownames(expression))) %>% pull(group_id)
        if (!is.null(start_cell)) {
            backbone_params$start.group.label <- groups_id %>% 
                filter(cell_id == start_cell) %>% pull(group_id)
        }
    }
    mst_tree <- do.call(cellTree::compute.backbone.tree, backbone_params)
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    cell_graph <- igraph::as_data_frame(mst_tree, "edges") %>% 
        dplyr::select(from, to, length = weight) %>% mutate(from = rownames(expression)[from], 
        to = rownames(expression)[to], directed = FALSE)
    to_keep <- igraph::V(mst_tree)$is.backbone %>% setNames(rownames(expression))
    tree <- cellTree:::.compute.tree.layout(mst_tree, ratio = 1)
    vertices <- igraph::as_data_frame(tree, "vertices") %>% as_data_frame()
    edges <- igraph::as_data_frame(tree, "edges") %>% as_data_frame()
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_cell_graph(cell_graph = cell_graph, to_keep = to_keep, 
            is_directed = FALSE, plot_vertices = vertices, plot_edges = edges) %>% 
        add_timings(timings = tl %>% add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')