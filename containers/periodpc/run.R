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
  expression,
  ndim = 3,
  maxit = 10
) {
  requireNamespace("stats")
  requireNamespace("princurve")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform PCA dimred
  dimred <- dimred(expression, method = "pca", ndim = ndim)

  # apply principal curve with periodic lowess smoother
  fit <- princurve::principal.curve(dimred, smoother = "periodic.lowess", maxit = maxit)

  # get pseudotime
  pseudotime <- fit$lambda %>% magrittr::set_names(rownames(expression))

  # construct segments
  path <- fit$s[fit$tag, , drop = FALSE]
  dimred_trajectory_segments <- cbind(
    path,
    path[c(seq(2, nrow(path)), 1), ,drop = F]
  ) %>%
    magrittr::set_colnames(c(paste0("from_", colnames(path)), paste0("to_", colnames(path))))

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_cyclic_trajectory(
    pseudotime = pseudotime
  ) %>% add_dimred(
    dimred = dimred,
    dimred_trajectory_segments = dimred_trajectory_segments
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")