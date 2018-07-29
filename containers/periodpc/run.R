library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)


#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_genes = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/periodpc/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression
#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

# perform PCA dimred
dimred <- dyndimred::dimred(expression, method = "pca", ndim = params$ndim)

# apply principal curve with periodic lowess smoother
fit <- princurve::principal_curve(dimred, smoother = "periodic.lowess", maxit = params$maxit)

# get pseudotime
pseudotime <- fit$lambda %>% magrittr::set_names(rownames(expression))

# construct segments
path <- fit$s[fit$ord, , drop = FALSE]
dimred_trajectory_segments <- cbind(
  path,
  path[c(seq(2, nrow(path)), 1), ,drop = F]
) %>%
  magrittr::set_colnames(c(paste0("from_", colnames(path)), paste0("to_", colnames(path))))

# TIMING: done with method
tl <- tl %>% add_timing_checkpoint("method_aftermethod")

# return output
model <- wrap_prediction_model(
  cell_ids = rownames(expression)
) %>% add_cyclic_trajectory(
  pseudotime = pseudotime
) %>% add_dimred(
  dimred = dimred,
  dimred_trajectory_segments = dimred_trajectory_segments
) %>% add_timings(
  timings = tl %>% add_timing_checkpoint("method_afterpostproc")
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")
