library(jsonlite)
library(readr)
library(dplyr)
library(purrr)


#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/identity/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts
dataset <- data$dataset

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# return output
output <- lst(
  cell_ids = dataset$cell_ids,
  milestone_ids = dataset$milestone_ids,
  milestone_network = dataset$milestone_network,
  divergence_regions = dataset$divergence_regions,
  progressions = dataset$progressions,
  timings = checkpoints
)


#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
