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
  expression,
  ndim = 2,
  dimred = "pca",
  component = 1
) {
  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  space <- dyndimred::dimred(expression, method = dimred, ndim = ndim)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotime = space[,component] %>% set_names(rownames(expression))
  ) %>% add_dimred(
    dimred = space
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")