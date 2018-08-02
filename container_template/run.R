library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(TEMPLATE)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/TEMPLATE/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

counts <- data$counts


# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))




# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())





# return output
model <- lst(
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
