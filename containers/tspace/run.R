library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(tSpace)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(model = "bifurcating", num_cells = 300) %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/tspace/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

expression <- data$expression

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

params <- list(
  K = 20,
  L = 15,
  D = "pearson_correlation",
  graph = 5,
  trajectories = 2,
  wp = 20,
  ground_truth = FALSE,
  weights = "exponential",
  dr = "pca"
)

fit <- tSpace::tSpace(
  expression,
  K = params$K,
  L = params$L,
  D = params$D,
  graph = params$graph,
  trajectories = params$trajectories,
  wp = params$wp,
  ground_truth = params$ground_truth,
  weights = params$weights,
  dr = params$dr
)


qplot(fit$ts_file$PC1, fit$ts_file$PC2) +
  geom_path(aes(fit$tspace_matrix[,1], fit$tspace_matrix[,2]))

fit$pca_embbeding$x

fit$tspace_matrix



# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())


output <- lst(
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
