#!/usr/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(Mpath)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/mpath/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts
groups_id <- data$groups_id

if (params$numcluster_null) {
  numcluster <- NULL
}

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# collect sample info
sample_info <- groups_id %>% select(cell_id, GroupID = group_id) %>% as.data.frame

# designate landmarks
landmark_cluster <- Mpath::landmark_designation(
  rpkmFile = t(counts),
  baseName = NULL,
  sampleFile = sample_info,
  distMethod = params$distMethod,
  method = params$method,
  numcluster = min(params$numcluster, nrow(counts) - 1),
  diversity_cut = params$diversity_cut,
  size_cut = params$size_cut,
  saveRes = FALSE
) %>%
  mutate_if(is.factor, as.character)

milestone_ids <- unique(landmark_cluster$landmark_cluster)

# catch situation where mpath only detects 1 landmark
if (length(milestone_ids) == 1) {
  warning("Mpath only detected one landmark")

  milestone_network <-
    data_frame(
      from = "M1",
      to = "M1",
      length = 0,
      directed = FALSE
    )
} else {
  # build network
  network <- Mpath::build_network(
    exprs = t(counts),
    baseName = NULL,
    landmark_cluster = landmark_cluster,
    distMethod = params$distMethod,
    writeRes = FALSE
  )

  # trim network
  trimmed_network <- Mpath::trim_net(
    nb12 = network,
    writeRes = FALSE
  )

  # create final milestone network
  class(trimmed_network) <- NULL
  milestone_network <- trimmed_network %>%
    reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
    mutate_if(is.factor, as.character) %>%
    filter(length > 0, from < to) %>%
    mutate(directed = FALSE)
}

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

grouping <-
  tibble::deframe(landmark_cluster)

output <- lst(
  cell_ids = names(grouping),
  grouping,
  milestone_network,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
