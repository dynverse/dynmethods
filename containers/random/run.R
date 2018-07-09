library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)



#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  counts,
  dummy_param = 0.5
) {
  num_milestones <- 15

  # generate network
  milestone_ids <- paste0("milestone_", seq_len(num_milestones))

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

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
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = cell_ids
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    divergence_regions = NULL
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
