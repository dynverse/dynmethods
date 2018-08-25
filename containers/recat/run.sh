#!/usr/local/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(reCAT)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/embeddr/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# run reCAT
result <- reCAT::bestEnsembleComplexTSP(
  test_exp = expression,
  TSPFold = params$TSPFold,
  beginNum = params$beginNum,
  endNum = params$endNum,
  base_cycle_range = seq(params$base_cycle_range_start, params$base_cycle_range_end),
  step_size = params$step_size,
  max_num = params$max_num,
  clustMethod = params$clustMethod,
  threads = 1,
  output = FALSE
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

pseudotime <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ] %>% set_names(rownames(expression))

# wrap
output <- lst(
  cell_ids = names(pseudotime),
  pseudotime,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
