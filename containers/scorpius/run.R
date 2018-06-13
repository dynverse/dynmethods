library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(SCORPIUS)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, ndim = 3, k = 4, distance_method = "spearman", 
    thresh = 0.001, maxit = 10, stretch = 0, smoother = "smooth.spline", 
    sparse = FALSE) 
{
    requireNamespace("SCORPIUS")
    if (k <= 1) {
        k <- NULL
    }
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    if (!sparse) {
        dist <- SCORPIUS::correlation_distance(expression, method = distance_method)
        space <- SCORPIUS::reduce_dimensionality(dist, ndim = ndim)
    }
    else {
        dist_fun <- function(x, y) SCORPIUS::correlation_distance(x, 
            y, method = distance_method)
        num_landmarks <- ifelse(nrow(expression) > 1000, 500, 
            nrow(expression))
        space <- SCORPIUS::reduce_dimensionality_landmarked(expression, 
            dist_fun = dist_fun, ndim = ndim, num_landmarks = num_landmarks)
    }
    traj <- SCORPIUS::infer_trajectory(space, k = k, thresh = thresh, 
        maxit = maxit, stretch = stretch, smoother = smoother)
    dimred_trajectory_segments <- cbind(traj$path[-nrow(traj$path), 
        , drop = FALSE] %>% magrittr::set_colnames(., paste0("from_comp_", 
        seq_along(colnames(.)))), traj$path[-1, , drop = FALSE] %>% 
        magrittr::set_colnames(., paste0("to_comp_", seq_along(colnames(.)))))
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_linear_trajectory(pseudotime = traj$time) %>% add_dimred(dimred = space, 
        dimred_trajectory_segments = dimred_trajectory_segments) %>% 
        add_timings(timings = tl %>% add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')