library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(ouija)
library(rstan)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_genes = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/embeddr/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# ouija assumes that a small number of marker genes is used, otherwise the method is too slow
expression <- expression[, data$features_id]

# write compiled instance of the stanmodel to HDD
rstan::rstan_options(auto_write = TRUE)

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# run ouija
oui <- ouija::ouija(
  x = expression,
  iter = params$iter,
  response_type = params$response_type,
  inference_type = params$inference_type,
  normalise_expression = params$normalise_expression
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# obtain the pseudotime
pseudotime <- ouija::map_pseudotime(oui) %>%
  setNames(rownames(expression))

# return output
output <- lst(
  pseudotime = pseudotime,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
