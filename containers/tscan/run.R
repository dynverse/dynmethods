library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(TSCAN)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "multifurcating") %>% c(., .$prior_information)
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

# order the cells
cds_order <- TSCAN::TSCANorder(cds_clus)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# process output
pseudotime <- set_names(seq_along(cds_order), cds_order)

dimred <- cds_clus$pcareducere

# return output
output <- lst(
  cell_ids = rownames(dimred),
  pseudotime,
  dimred,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
