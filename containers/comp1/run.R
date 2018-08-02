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
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "binary_tree") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/comp1/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

expression <- data$expression

# TIMING: done with preproc
tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

space <- dyndimred::dimred(expression, method = params$dimred, ndim = params$ndim)

# TIMING: done with method
tl <- tl %>% add_timing_checkpoint("method_aftermethod")

# return output
model <- wrap_prediction_model(
  cell_ids = rownames(expression)
) %>% add_linear_trajectory(
  pseudotime = space[,params$component] %>% set_names(rownames(expression))
) %>% add_dimred(
  dimred = space
) %>% add_timings(
  timings = tl %>% add_timing_checkpoint("method_afterpostproc")
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
