library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

source("SupplementaryMethods/Waterfall.R")


#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' source("~/Downloads/SupplementaryMethods/Waterfall.R")
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/waterfall/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# run waterfall
ps <- pseudotimeprog.foo(t(expression), k = params$num_clusters, color = rep("black", nrow(expression)))

dimred <- ps[,colnames(ps) != "pseudotime", drop = FALSE]

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# return output
output <- lst(
  pseudotime = set_names(ps$pseudotime, rownames(ps)),
  dimred = as.matrix(dimred),
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
