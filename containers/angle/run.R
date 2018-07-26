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
#' data <- dyntoy::generate_dataset(model = "cyclic") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/angle/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

expression <- data$expression


# TIMING: done with preproc
tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

# perform PCA dimred
space <- dyndimred::dimred(expression, method = params$dimred, ndim = 2)

# transform to pseudotime using atan2
pseudotime <- atan2(space[,2], space[,1]) / 2 / pi + .5

# TIMING: done with method
tl <- tl %>% add_timing_checkpoint("method_aftermethod")

# return output
model <- wrap_prediction_model(
  cell_ids = rownames(expression)
) %>% add_cyclic_trajectory(
  pseudotime = pseudotime,
  do_scale_minmax = FALSE
) %>% add_dimred(
  dimred = space
) %>% add_timings(
  timings = tl %>% add_timing_checkpoint("method_afterpostproc")
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
