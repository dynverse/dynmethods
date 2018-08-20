library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(monocle)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/monocle_ddrtree/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts
#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# just in case
if (is.factor(params$norm_method)) {
  params$norm_method <- as.character(params$norm_method)
}

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# load in the new dataset
pd <- Biobase::AnnotatedDataFrame(data.frame(row.names = rownames(counts)))
fd <- Biobase::AnnotatedDataFrame(data.frame(row.names = colnames(counts), gene_short_name = colnames(counts)))
cds <- monocle::newCellDataSet(t(counts), pd, fd)

# estimate size factors and dispersions
cds <- BiocGenerics::estimateSizeFactors(cds)
cds <- BiocGenerics::estimateDispersions(cds)

# if low # cells or features -> https://github.com/cole-trapnell-lab/monocle-release/issues/26
# this avoids the error "initial centers are not distinct."
if (ncol(counts) < 500 || nrow(counts) < 500) {
  params$auto_param_selection <- FALSE
}

# reduce the dimensionality
cds <- monocle::reduceDimension(
  cds,
  max_components = params$max_components,
  reduction_method = params$reduction_method,
  norm_method = params$norm_method,
  auto_param_selection = params$auto_param_selection
)

# order the cells
cds <- monocle::orderCells(cds)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# extract the igraph and which cells are on the trajectory
gr <- cds@auxOrderingData[[params$reduction_method]]$pr_graph_cell_proj_tree
to_keep <- setNames(rep(TRUE, nrow(counts)), rownames(counts))

# convert to milestone representation
cell_graph <- igraph::as_data_frame(gr, "edges") %>% mutate(directed = FALSE)

if ("weight" %in% colnames(cell_graph)) {
  cell_graph <- cell_graph %>% rename(length = weight)
} else {
  cell_graph <- cell_graph %>% mutate(length = 1)
}

cell_graph <- cell_graph %>% select(from, to, length, directed)

dimred <- t(cds@reducedDimS)
colnames(dimred) <- paste0("Comp", seq_len(ncol(dimred)))

# wrap output
output <- lst(
  cell_ids = rownames(dimred),
  cell_graph,
  to_keep,
  dimred,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
