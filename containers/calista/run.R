library(dplyr)
library(purrr)
library(readr)
library(feather)

# calista NEEDS to be in the CALISTA-R folder while at the same time
# requiring to be able to write files there
file.copy("/CALISTA/CALISTA-R", "/ti/workspace", recursive = TRUE)
setwd("/ti/workspace/CALISTA-R")
source("R/initialization.R")

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(model = "cyclic") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/calista/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)
#' setwd("~/Downloads/CALISTA/CALISTA-R/")

expression <- data$expression
file_loc <- tempfile(pattern = "expression.csv")

data_df <- data.frame(
  row.names = NULL,
  expression,
  check.names = FALSE
)
write_csv(data_df, file_loc)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# Prepare CALISTA for work
INPUTS <- list()
# Specify data types and settings for pre-processing
INPUTS$data_location=file_loc
INPUTS$data_type=1; # Single-cell RT-qPCR CT data
INPUTS$format_data=3; # Rows= cells and Columns= genes (no time/stage info)
INPUTS$data_selection= integer(); # Include data from all time points
INPUTS$perczeros_genes=100; # Remove genes with > 100# of zeros
INPUTS$perczeros_cells=100; # Remove cells with 100# of zeros
INPUTS$cells_2_cut=0; # No manual removal of cells
INPUTS$perc_top_genes=100; # Retain only top X the most variable genes with X=min(200, INPUTS$perc_top_genes * num of cells/100, num of genes)
# Specify single-cell clustering settings
INPUTS$optimize=1; # The number of cluster is known a priori
INPUTS$parallel=0; # Use parallelization option
INPUTS$cluster='kmedoids'; # Use k-medoids in consensus clustering
INPUTS$Cluster='kmedoids'; # Use k-medoids in consensus clustering
# Specify transition genes settings
INPUTS$thr_transition_genes=50; # Set threshold for transition genes determination to

# put parameters into INPUTS
INPUTS$runs=params$runs; # Perform 50 independent runs of greedy algorithm
INPUTS$max_iter=params$max_iter; # Limit the number of iterations in greedy algorithm to 100

# % Upload and pre-process data
DATA=import_data(INPUTS)

checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# %% *** 2-SINGLE-CELL CLUSTERING ***
# %
# % Please check comments in 'CALISTA_clustering_main' for more information.
# %

Results=list()
CALISTA_clustering_main_results=CALISTA_clustering_main(DATA,INPUTS)
Results=CALISTA_clustering_main_results$Results
DATA=CALISTA_clustering_main_results$DATA
INPUTS=CALISTA_clustering_main_results$INPUTS

###cluster cutting
cluster_cut=0L

# %% *** 3-RECONSTRUCTION OF LINEAGE PROGRESSION ***
# %
# % Please check comments in 'CALISTA_transition_main' for more information.
Results=CALISTA_transition_main(DATA,Results)

# %% *** 4-DETERMINATION OF TRANSITION GENES ***
# %
# % Please check comments in 'CALISTA_transition_genes_main' for more information.
Results=CALISTA_transition_genes_main(DATA,INPUTS,Results)

# %% *** 5-PSEUDOTEMPORAL ORDERING OF CELLS ***
# %
# % Please check comments in 'CALISTA_ordering_main' for more information.
#
Results=CALISTA_ordering_main(DATA,Results)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# CREATE OUTPUT
milestone_network <- inner_join(data.frame(Results$TRANSITION$nodes_connection), data.frame(Results$cluster_distance), by = c("X1", "X2")) %>%
  dplyr::select(from = X1, to = X2, length = X3) %>%
  mutate(directed = FALSE)

cluster_probs <- Results$clustering_struct$all$all$clusterprobabilities

progressions <- map_df(seq_len(nrow(cluster_probs)), function(cell_i) {
  if (!any(is.finite(cluster_probs[cell_i, , drop = TRUE]))) {
    NULL
  } else {
    milestone_network %>%
      mutate(
        from_prob = cluster_probs[cell_i, from, drop = TRUE],
        to_prob = cluster_probs[cell_i, to, drop = TRUE],
        sum = from_prob + to_prob,
        percentage = to_prob / sum,
        cell_id = rownames(expression)[[cell_i]]
      ) %>%
      arrange(desc(sum)) %>%
      slice(1) %>%
      dplyr::select(cell_id, from, to, percentage)
  }
})

milestone_network <- milestone_network %>%
  mutate(from = paste0("M", from), to = paste0("M", to))

progressions <- progressions %>%
  mutate(from = paste0("M", from), to = paste0("M", to))

output <- lst(
  cell_ids = unique(progressions$cell_id),
  progressions,
  milestone_network,
  timings = checkpoints
)

write_rds(output, "/ti/output/output.rds")
