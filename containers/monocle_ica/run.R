library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(BiocGenerics)
library(igraph)
library(Biobase)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (counts, groups_n = NULL, reduction_method = "ICA", 
    max_components = 2, norm_method = "vstExprs", auto_param_selection = TRUE) 
{
    requireNamespace("monocle")
    requireNamespace("BiocGenerics")
    requireNamespace("igraph")
    requireNamespace("Biobase")
    requireNamespace("Matrix")
    if (is.factor(norm_method)) 
        norm_method <- as.character(norm_method)
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    pd <- Biobase::AnnotatedDataFrame(data.frame(row.names = rownames(counts)))
    fd <- Biobase::AnnotatedDataFrame(data.frame(row.names = colnames(counts), 
        gene_short_name = colnames(counts)))
    cds <- monocle::newCellDataSet(t(counts), pd, fd)
    cds <- BiocGenerics::estimateSizeFactors(cds)
    cds <- BiocGenerics::estimateDispersions(cds)
    cds <- monocle::reduceDimension(cds, max_components = max_components, 
        reduction_method = reduction_method, norm_method = norm_method, 
        auto_param_selection = auto_param_selection)
    cds <- monocle::orderCells(cds, num_paths = groups_n)
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    gr <- if (reduction_method == "DDRTree") {
        cds@auxOrderingData[[reduction_method]]$pr_graph_cell_proj_tree
    }
    else if (reduction_method == "ICA") {
        cds@auxOrderingData[[reduction_method]]$cell_ordering_tree
    }
    to_keep <- setNames(rep(TRUE, nrow(counts)), rownames(counts))
    cell_graph <- igraph::as_data_frame(gr, "edges") %>% mutate(directed = FALSE)
    if ("weight" %in% colnames(cell_graph)) {
        cell_graph <- cell_graph %>% rename(length = weight)
    }
    else {
        cell_graph <- cell_graph %>% mutate(length = 1)
    }
    cell_graph <- cell_graph %>% select(from, to, length, directed)
    if (exists("postprocess_monocle_cds")) {
        plot_data <- postprocess_monocle_cds(cds)
    }
    else {
        plot_data <- NULL
    }
    wrap_prediction_model(cell_ids = rownames(counts)) %>% add_cell_graph(cell_graph = cell_graph, 
        to_keep = to_keep, plot_data = plot_data, reduction_method = reduction_method) %>% 
        add_timings(timings = tl %>% add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')