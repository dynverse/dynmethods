#!/usr/local/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(TSCAN)
library(igraph)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_features = 300, model = "multifurcating") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/tscan/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# process clusternum
clusternum <- seq(params$clusternum_lower, params$clusternum_upper, 1)

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# preprocess counts
cds_prep <- TSCAN::preprocess(
  t(as.matrix(counts)),
  takelog = TRUE,
  logbase = 2,
  pseudocount = 1,
  clusternum = NULL,
  minexpr_value = params$minexpr_value,
  minexpr_percent = params$minexpr_percent,
  cvcutoff = params$cvcutoff
)

# cluster the data
cds_clus <- TSCAN::exprmclust(
  cds_prep,
  clusternum = clusternum,
  modelNames = params$modelNames,
  reduce = TRUE
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# process output
milestone_network <- cds_clus$MSTtree %>%
  igraph::as_data_frame() %>%
  rename(length = weight) %>%
  mutate(directed = FALSE)
dimred <- cds_clus$pcareduceres
dimred_milestones <- cds_clus$clucenter
rownames(dimred_milestones) <- as.character(seq_len(nrow(dimred_milestones)))
colnames(dimred_milestones) <- colnames(dimred)

# return output
output <- lst(
  cell_ids = rownames(dimred),
  milestone_ids = rownames(dimred_milestones),
  milestone_network = milestone_network,
  dimred_milestones,
  dimred,
  grouping = cds_clus$clusterid,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
