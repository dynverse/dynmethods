library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(SCORPIUS)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/scorpius/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# use k <= 1 to turn off clustering
if (params$k <= 1) {
  params$k <- NULL
}

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

space <- SCORPIUS::reduce_dimensionality(
  x = expression,
  dist_fun = function(x, y = NULL) dynutils::calculate_distance(x = x, y = y, method = params$distance_method),
  landmark_method = ifelse(params$sparse, "naive", "none"),
  ndim = params$ndim,
  num_landmarks = ifelse(nrow(expression) > 1000, 500, nrow(expression))
)

# infer a trajectory through the data
traj <- SCORPIUS::infer_trajectory(
  space,
  k = params$k,
  thresh = params$thresh,
  maxit = params$maxit,
  stretch = params$stretch,
  smoother = params$smoother
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# convert trajectory to segments
dimred_trajectory_segments <-
  cbind(
    traj$path[-nrow(traj$path), , drop = FALSE] %>%
      magrittr::set_colnames(., paste0("from_comp_", seq_len(ncol(.)))),
    traj$path[-1, , drop = FALSE] %>%
      magrittr::set_colnames(., paste0("to_comp_", seq_len(ncol(.))))
  )

# return output
output <- lst(
  pseudotime = traj$time,
  dimred = space,
  dimred_trajectory_segments = dimred_trajectory_segments,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
