library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(TSCAN)
library(igraph)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

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
cluster_network <- cds_clus$MSTtree %>%
  igraph::as_data_frame() %>%
  rename(length = weight) %>%
  mutate(directed = FALSE)
sample_space <- cds_clus$pcareduceres
cluster_space <- cds_clus$clucenter
rownames(cluster_space) <- as.character(seq_len(nrow(cluster_space)))
colnames(cluster_space) <- colnames(sample_space)

# return output
output <- lst(
  milestone_ids = rownames(cluster_space),
  milestone_network = cluster_network,
  dimred_milestones = cluster_space,
  dimred = sample_space,
  milestone_assignment_cells = cds_clus$clusterid,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
