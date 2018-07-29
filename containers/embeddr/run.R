library(jsonlite)
library(readr)
library(dplyr)
library(purrr)
requireNamespace("scaterlegacy")
requireNamespace("embeddr")

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_genes = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/embeddr/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

counts <- data$counts

# calculate nn param
nn <- max(round(log(nrow(counts)) * params$nn_pct), 9)

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# load data in scaterlegacy
sce <- scaterlegacy::newSCESet(countData = t(counts))

# run embeddr
sce <- embeddr::embeddr(
  sce,
  kernel = params$kernel,
  metric = params$metric,
  nn = nn,
  eps = params$eps,
  t = params$t,
  symmetrize = params$symmetrize,
  measure_type = params$measure_type,
  p = params$ndim
)

# fit pseudotime
sce <- embeddr::fit_pseudotime(
  sce,
  thresh = params$thresh,
  maxit = params$maxit,
  stretch = params$stretch,
  smoother = params$smoother
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# construct milestone network
pseudotime <- as(sce@phenoData, "data.frame")$pseudotime %>%
  setNames(rownames(counts))

# creating extra output for visualisation purposes
dimred_cells <- sce@reducedDimension

traj <- as(sce@phenoData, "data.frame") %>%
  dplyr::arrange(pseudotime) %>%
  select(starts_with("trajectory_")) %>%
  as.matrix()

dimred_trajectory_segments <- cbind(
  traj[-nrow(traj), , drop = F],
  traj[-1, , drop = F]
)
colnames(dimred_trajectory_segments) <- c(
  paste0("from_comp_", seq_len(ncol(dimred_cells))),
  paste0("to_comp_", seq_len(ncol(dimred_cells)))
)

# return output
model <- lst(
  pseudotime = pseudotime,
  dimred = dimred_cells,
  dimred_trajectory_segments = dimred_trajectory_segments,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
