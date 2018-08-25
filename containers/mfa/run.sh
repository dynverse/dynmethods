#!/usr/local/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

library(mfa)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 100, num_features = 101, model = "bifurcating") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/mfa/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# perform MFA
m <- mfa::mfa(
  y = expression,
  b = params$b,
  iter = params$iter,
  thin = params$thin,
  zero_inflation = params$zero_inflation,
  pc_initialise = params$pc_initialise,
  prop_collapse = params$prop_collapse,
  scale_input = params$scale_input
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# obtain results
ms <- summary(m) %>%
  mutate(cell_id = rownames(expression)) %>%
  select(cell_id, everything()) %>%
  group_by(cell_id) %>%
  mutate(
    branch = paste0("M", branch),
    branch_certainty = branch_certainty / sum(branch_certainty)
  ) %>%
  ungroup()

end_state_probabilities <- ms %>%
  select(cell_id, branch, branch_certainty) %>%
  spread(branch, branch_certainty, fill = 0)

pseudotime <-
  ms %>%
  group_by(cell_id) %>%
  summarise(pseudotime = sum(branch_certainty * pseudotime)) %>%
  deframe()


# return output
output <- lst(
  cell_ids = names(pseudotime),
  end_state_probabilities,
  pseudotime,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
