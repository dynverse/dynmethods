library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(phenopath)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_genes = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/embeddr/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# run phenopath
fit <- phenopath::phenopath(
  exprs_obj = expression,
  x = rep(1, nrow(expression)),
  elbo_tol = 1e-6,
  thin = params$thin,
  z_init = ifelse(params$z_init == "random", "random", as.numeric(params$z_init)),
  model_mu = params$model_mu,
  scale_y = params$scale_y
)
pseudotime <- phenopath::trajectory(fit) %>%
  setNames(rownames(expression))

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# return output
output <- lst(
  pseudotime,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
