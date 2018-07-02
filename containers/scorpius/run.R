library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(SCORPIUS)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  expression,
  ndim = 3,
  k = 4,
  distance_method = "spearman",
  thresh = 0.001,
  maxit = 10,
  stretch = 0,
  smoother = "smooth_spline",
  sparse = FALSE
) {
  requireNamespace("SCORPIUS")

  # if k is too low, turn off clustering
  if (k <= 1) {
    k <- NULL
  }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  space <- SCORPIUS::reduce_dimensionality(
    x = expression,
    dist_fun = function(x, y = NULL) dynutils::calculate_distance(x = x, y = y, method = distance_method),
    landmark_method = ifelse(sparse, "naive", "none"),
    ndim = ndim,
    num_landmarks = ifelse(nrow(expression) > 1000, 500, nrow(expression))
  )

  # infer a trajectory through the data
  traj <- SCORPIUS::infer_trajectory(
    space,
    k = k,
    thresh = thresh,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother
  )

  # convert trajectory to segments
  dimred_trajectory_segments <-
    cbind(
      traj$path[-nrow(traj$path), , drop = FALSE] %>% magrittr::set_colnames(., paste0("from_comp_", seq_along(colnames(.)))),
      traj$path[-1, , drop = FALSE] %>% magrittr::set_colnames(., paste0("to_comp_", seq_along(colnames(.))))
    )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotime = traj$time
  ) %>% add_dimred(
    dimred = space,
    dimred_trajectory_segments = dimred_trajectory_segments
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
