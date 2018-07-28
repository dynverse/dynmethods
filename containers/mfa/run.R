library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(mfa)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_genes = 300, model = "linear") %>% c(., .$prior_information)
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
  select(cell_id, everything())

# create milestone network
milestone_network <- data_frame(
  from = "M0",
  to = paste0("M", seq_len(params$b)),
  length = 1,
  directed = TRUE
)

# create progressions
progressions <- with(ms, data_frame(
  cell_id = rownames(expression),
  from = "M0",
  to = paste0("M", branch),
  percentage = (pseudotime - min(pseudotime)) / (max(pseudotime) - min(pseudotime))
))

# return output
output <- lst(
  milestone_network,
  progressions,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
