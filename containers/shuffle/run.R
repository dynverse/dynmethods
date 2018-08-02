library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/embeddr/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts
dataset <- data$dataset

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# permute cell labels
allcells <- rownames(counts)
mapper <- setNames(sample(allcells), allcells)
progressions <- dataset$progressions %>% mutate(
  cell_id = mapper[cell_id]
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

output <- lst(
  milestone_network = dataset$milestone_network,
  progressions = progressions,
  divergence_regions = dataset$divergence_regions,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
