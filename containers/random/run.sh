#!/usr/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

counts <- data$counts

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

num_milestones <- 15

# generate network
milestone_ids <- paste0("milestone_", seq_len(num_milestones))

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

gr <- igraph::ba.game(num_milestones)
milestone_network <- igraph::as_data_frame(gr) %>%
  mutate(
    from = paste0("milestone_", from),
    to = paste0("milestone_", to),
    length = 1,
    directed = FALSE
  )

# put cells on random edges of network
cell_ids <- rownames(counts)

progressions <- data.frame(
  cell_id = cell_ids,
  milestone_network[sample.int(nrow(milestone_network), length(cell_ids), replace = TRUE), 1:2],
  percentage = runif(length(cell_ids)),
  stringsAsFactors = FALSE
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# return output
output <- lst(
  cell_ids = cell_ids,
  milestone_ids = milestone_ids,
  milestone_network = milestone_network,
  progressions = progressions,
  divergence_regions = NULL,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
