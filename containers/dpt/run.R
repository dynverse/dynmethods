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
#' data <- dyntoy::generate_dataset(model = "tree") %>% c(., .$prior_information)
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

cell_ids <- rownames(expression)

# construct grouping
grouping <- dpt@branch[,1] %>%
  ifelse(is.na(.), 0, .) %>%
  as.character() %>%
  paste0("Tip", .) %>%
  set_names(cell_ids)

group_ids <- sort(unique(grouping))

# collect distances from tips
tip_dists <- dpt[,tips] %>%
  magrittr::set_colnames(., paste0("Tip", seq_len(ncol(.)))) %>%
  magrittr::set_rownames(cell_ids)

# calculate progressions
outs <- map(
  group_ids,
  function(gid) {
    cat("Processing ", gid, "\n", sep = "")
    cixs <- which(grouping == gid)
    cids <- cell_ids[cixs]
    if (length(cids) > 0) {
      if (gid == "Tip0") {
        progr <- data_frame(
          cell_id = cids,
          from = gid,
          to = sample(setdiff(group_ids, gid), length(cids), replace = TRUE),
          percentage = 0
        )
        list(progr = progr, milnet = NULL)
      } else {
        # calculate min dist of gid to all other cells
        max_range <- min(tip_dists[-cixs, gid])

        # calculate percentage value of cells in gid
        percentage <- 1 - pmin(tip_dists[cixs, gid] / max_range, 1)
        progr <- data_frame(cell_id = cids, from = "Tip0", to = gid, percentage = percentage)
        milnet <- data_frame(from = "Tip0", to = gid, length = max_range, directed = FALSE)
        list(progr = progr, milnet = milnet)
      }
    } else {
      list()
    }
  }
)
progressions <- map_df(outs, ~ .$progr)
milestone_network <- map_df(outs, ~ .$milnet)

# collect dimred
dimred <- dm@eigenvectors %>%
  magrittr::set_colnames(., paste0("Comp", seq_len(ncol(.)))) %>%
  magrittr::set_rownames(cell_ids)

# return output
output <- lst(
  cell_ids,
  milestone_ids = group_ids,
  milestone_network,
  progressions,
  group_ids,
  grouping,
  dimred,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
