library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)



#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  counts,
  dataset,
  dummy_param = 0.5
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = dataset$cell_ids,
    cell_info = dataset$cell_info
  ) %>% add_trajectory(
    milestone_ids = dataset$milestone_ids,
    milestone_network = dataset$milestone_network,
    divergence_regions = dataset$divergence_regions,
    progressions = dataset$progressions
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")