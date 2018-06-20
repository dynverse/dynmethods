library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library()

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  counts,
  dummy_param = 0.5
) {
  stop("This control method always errors.")
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')