#!/usr/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(destiny)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(model = "bifurcating") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/dpt/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

expression <- data$expression

start_cell <-
  if (!is.null(data$start_id)) {
    sample(data$start_id, 1)
  } else {
    NULL
  }

# create n_local vector
n_local <- seq(params$n_local_lower, params$n_local_upper, by = 1)

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# run diffusion maps
dm <- destiny::DiffusionMap(
  data = expression,
  sigma = params$sigma,
  distance = params$distance,
  n_eigs = params$ndim,
  density_norm = params$density_norm,
  n_local = n_local,
  vars = params$features_id
)

# run DPT
if (!is.null(data$start_cell)) {
  tips <- which(rownames(expression) %in% start_cell)
} else {
  tips <- destiny::random_root(dm)
}
dpt <- destiny::DPT(
  dm,
  w_width = params$w_width,
  tips = tips
)

# find DPT tips
tips <- destiny::tips(dpt)
tip_names <- rownames(expression)[tips]

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# retrieve dimred
dimred_cells <- dpt@dm@eigenvectors %>% magrittr::set_rownames(rownames(expression)) %>% as.matrix

# get cluster assignment
grouping <- dpt@branch[,1] %>%
  ifelse(is.na(.), 0, .) %>%
  as.character()
branches <- sort(unique(grouping))

# calculate cluster medians
dimred_milestones <- t(sapply(branches, function(br) colMeans(dimred_cells[grouping == br,,drop=F])))

# create star network
milestone_network <- data_frame(
  from = "0",
  to = setdiff(branches, "0"),
  length = sqrt(rowMeans((dimred_milestones[from,] - dimred_milestones[to,])^2)),
  directed = FALSE
)

output <- lst(
  cell_ids = rownames(dimred_cells),
  milestone_network,
  dimred_milestones,
  dimred = dimred_cells,
  grouping = grouping,
  tips,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
