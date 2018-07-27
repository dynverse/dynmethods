library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(destiny)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(model = "bifurcating") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/kbranches/definition.yml")$parameters %>%
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

#keep the first 2 diffusion components
input_dat <- destiny::as.data.frame(dm)[,1:2]

# Clustingering probably doesn't lead to information that I can enter into DPT
# clust <- kbranches::kbranches.global(input_dat, Kappa = 3)

Dmat <- kbranches::compute_all_distances(input_dat)
res <- kbranches::kbranches.local(input_dat = input_dat, Dmat = Dmat)
tips <- kbranches::identify_regions(
  input_dat = input_dat,
  gap_scores = res$gap_scores,
  Dist = Dmat,
  smoothing_region = 5,
  smoothing_region_thresh = 5,
  mode = 'tip'
)
tips <- map_int(seq_len(tips$num_clusters), ~ sample(which(tips$cluster == .), 1))

dpt <- destiny::DPT(dm, tips = tips)


# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

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

output <- lst(
  milestone_network,
  dimred_milestones,
  dimred = dimred_cells,
  milestone_assignment_cells,
  tips,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
