library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(phenopath)
library(dyndimred)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  expression,
  thin = 40,
  z_init = "1",
  model_mu = FALSE,
  scale_y = TRUE
) {
  requireNamespace("phenopath")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run phenopath
  fit <- phenopath::phenopath(
    exprs_obj = expression,
    x = rep(1, nrow(expression)),
    elbo_tol = 1e-6,
    thin = thin,
    z_init = ifelse(z_init == "random", "random", as.numeric(z_init)),
    model_mu = model_mu,
    scale_y = scale_y
  )
  pseudotime <- phenopath::trajectory(fit) %>%
    setNames(rownames(expression))

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # run pca for visualisation purposes
  space <- dyndimred::dimred(expression, method = "pca", ndim = 2)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotime = pseudotime
  ) %>% add_dimred(
    dimred = space
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")