library(jsonlite)
library(readr)
library(dplyr)
library(purrr)
library(dyndimred)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(model = "cyclic") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/angle/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

expression <- data$expression

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# perform PCA dimred
space <- dyndimred::dimred(expression, method = params$dimred, ndim = 2)

# transform to pseudotime using atan2
pseudotime <- atan2(space[,2], space[,1]) / 2 / pi + .5

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# return output
output <- lst(
  cell_ids = rownames(expression),
  pseudotime = pseudotime,
  do_scale_minmax = FALSE,
  dimred = space,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
