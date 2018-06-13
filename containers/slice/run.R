library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(SLICE)
library(igraph)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, groups_id = NULL, features_id = NULL, lm.method = "clustering", 
    model.type = "tree", ss.method = "all", ss.threshold = 0.25, 
    community.method = "louvain", cluster.method = "kmeans", 
    k = 0, k.max = 10, B = 100, k.opt.method = "firstmax") 
{
    requireNamespace("SLICE")
    requireNamespace("igraph")
    if (k == 0) {
        k <- NULL
    }
    if (!is.null(groups_id)) {
        cellidentity <- groups_id %>% slice(match(rownames(expression), 
            cell_id)) %>% pull(group_id) %>% factor()
    }
    else {
        cellidentity <- factor(rep(1, nrow(expression)))
    }
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    sc <- SLICE::construct(exprmatrix = as.data.frame(t(expression)), 
        cellidentity = cellidentity)
    num_genes <- ncol(expression)
    km <- matrix(runif(num_genes * num_genes), ncol = num_genes, 
        dimnames = list(colnames(expression), colnames(expression)))
    sc <- SLICE::getEntropy(sc, km = km)
    sc <- SLICE::getRDS(sc, method = "pca", num_dim = 2, log.base = 2, 
        do.center = TRUE, do.scale = FALSE, use.cor = TRUE, min.var = 0, 
        min.cells = 0, genes.use = features_id)
    sc <- SLICE::getLineageModel(sc, model.type = model.type, 
        ss.method = ss.method, ss.threshold = ss.threshold, community.method = community.method, 
        cluster.method = cluster.method, k = k, k.max = k.max, 
        B = B, k.opt.method = k.opt.method, do.plot = FALSE)
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    lin_model <- sc@model$lineageModel
    milestone_ids <- names(igraph::V(lin_model))
    milestone_network <- lin_model %>% igraph::as_data_frame() %>% 
        rename(length = weight) %>% mutate(directed = TRUE)
    pseudotime <- map_df(seq_len(nrow(milestone_network)), function(i) {
        from <- milestone_network[i, 1]
        to <- milestone_network[i, 2]
        sc_tmp <- SLICE::getTrajectories(sc, method = "pc", start = match(from, 
            milestone_ids), end = match(to, milestone_ids), do.plot = FALSE, 
            do.trim = FALSE)
        sc_tmp@transitions[[1]]$i.pseudotime %>% tibble::rownames_to_column("cell_id") %>% 
            mutate(from = from, to = to) %>% select(cell_id, 
            from, to, percentage = ptime)
    })
    progressions <- sc@model$cells.df %>% tibble::rownames_to_column("cell_id") %>% 
        slice(match(rownames(expression), cell_id)) %>% mutate(state = paste0("slice.ss.", 
        slice.state)) %>% select(cell_id, state) %>% right_join(pseudotime, 
        by = "cell_id") %>% filter((state == from) | (state == 
        to)) %>% group_by(cell_id) %>% arrange(percentage) %>% 
        slice(1) %>% select(-state) %>% ungroup()
    cells.df <- sc@model$cells.df
    edge.df <- igraph::get.edgelist(lin_model) %>% as.data.frame() %>% 
        mutate(ix = match(V1, rownames(cells.df)), iy = match(V2, 
            rownames(cells.df)), src.x = cells.df$x[ix], src.y = cells.df$y[ix], 
            dst.x = cells.df$x[iy], dst.y = cells.df$y[iy])
    dimred <- cells.df %>% tibble::rownames_to_column("cell_id") %>% 
        filter(slice.realcell == 1) %>% tibble::column_to_rownames("cell_id") %>% 
        select(x, y) %>% as.matrix()
    dimred_milestones <- cells.df %>% tibble::rownames_to_column("cell_id") %>% 
        filter(slice.realcell == 0) %>% tibble::column_to_rownames("cell_id") %>% 
        select(x, y) %>% as.matrix()
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_trajectory(milestone_ids = milestone_ids, milestone_network = milestone_network, 
            progressions = progressions, divergence_regions = NULL) %>% 
        add_dimred(dimred = dimred, dimred_milestones = dimred_milestones, 
            dimred_trajectory_segments = edge.df[, c("src.x", 
                "src.y", "dst.x", "dst.y")] %>% mutate_all(as.numeric) %>% 
                as.matrix %>% magrittr::set_colnames(c("from_comp_1", 
                "from_comp_2", "to_comp_1", "to_comp_2")), cells.df = cells.df, 
            edge.df = edge.df) %>% add_timings(timings = tl %>% 
        add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')