library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(mfa)
library(dyndimred)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  expression,
  b = 2,
  iter = 2000,
  thin = 1,
  zero_inflation = FALSE,
  pc_initialise = 1,
  prop_collapse = 0,
  scale_input = TRUE
) {
  requireNamespace("mfa")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform MFA
  m <- mfa::mfa(
    y = expression,
    b = b,
    iter = iter,
    thin = thin,
    zero_inflation = zero_inflation,
    pc_initialise = pc_initialise,
    prop_collapse = prop_collapse,
    scale_input = scale_input
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # obtain results
  ms <- summary(m) %>%
    mutate(cell_id = rownames(expression)) %>%
    select(cell_id, everything())

  # create milestone network
  milestone_ids <- paste0("M", seq(0, b))
  milestone_network <- data_frame(
    from = "M0",
    to = paste0("M", seq_len(b)),
    length = 1,
    directed = TRUE
  )

  # create progressions
  progressions <- with(ms, data_frame(
    cell_id = rownames(expression),
    from = "M0",
    to = paste0("M", branch),
    percentage = dynutils::scale_minmax(pseudotime)
  ))

  # pca for visualisation only
  pca_out <- dyndimred::dimred(expression, method = "pca", ndim = 2)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression),
    cell_info = ms
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    divergence_regions = NULL
  ) %>% add_dimred(
    dimred = pca_out
  ) %>% add_timings(
    tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')