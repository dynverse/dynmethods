library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(dynutils)
library(reshape2)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  expression,
  start_id = NULL,
  features_id = NULL,
  sigma = "local",
  distance = "euclidean",
  ndim = 20,
  density_norm = TRUE,
  n_local_lower = 5,
  n_local_upper = 7,
  w_width = 0.1
) {
  requireNamespace("destiny")

  start_cell <-
    if (!is.null(start_id)) {
      sample(start_id, 1)
    } else {
      NULL
    }

  # create n_local vector
  n_local <- seq(n_local_lower, n_local_upper, by = 1)

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run diffusion maps
  dm <- destiny::DiffusionMap(
    data = expression,
    sigma = sigma,
    distance = distance,
    n_eigs = ndim,
    density_norm = density_norm,
    n_local = n_local,
    vars = features_id
  )

  # run DPT
  dpt_params <- lst(dm, w_width)
  if (!is.null(start_cell)) {
    dpt_params$tips <- which(rownames(expression) %in% start_cell)
  }
  dpt <- do.call(destiny::DPT, dpt_params)

  # find DPT tips
  tips <- destiny::tips(dpt)
  tip_names <- rownames(expression)[tips]

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # retrieve dimred
  dimred_cells <- dpt@dm@eigenvectors %>% magrittr::set_rownames(rownames(expression)) %>% as.matrix

  # get cluster assignment
  milestone_assignment_cells <- dpt@branch[,1] %>%
    ifelse(is.na(.), 0, .) %>%
    as.character()
  branches <- sort(unique(milestone_assignment_cells))

  # calculate cluster medians
  dimred_milestones <- t(sapply(branches, function(br) colMeans(dimred_cells[milestone_assignment_cells == br,,drop=F])))

  # create star network
  milestone_network <- data_frame(
    from = "0",
    to = setdiff(branches, "0"),
    length = sqrt(rowMeans((dimred_milestones[from,] - dimred_milestones[to,])^2)),
    directed = TRUE
  )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_dimred_projection(
    milestone_ids = rownames(dimred_milestones),
    milestone_network = milestone_network,
    dimred_milestones = dimred_milestones,
    dimred = dimred_cells,
    milestone_assignment_cells = milestone_assignment_cells,
    tips = tip_names
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")