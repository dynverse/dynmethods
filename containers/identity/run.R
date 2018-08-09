library(dynwrap)
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
#' params <- yaml::read_yaml("containers/identity/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts
dataset <- data$dataset

# TIMING: done with preproc
tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

# TIMING: done with method
tl <- tl %>% add_timing_checkpoint("method_aftermethod")

# return output
model <-
  wrap_prediction_model(
    cell_ids = dataset$cell_ids
  ) %>% add_trajectory(
    milestone_ids = dataset$milestone_ids,
    milestone_network = dataset$milestone_network,
    divergence_regions = dataset$divergence_regions,
    progressions = dataset$progressions
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )


#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
