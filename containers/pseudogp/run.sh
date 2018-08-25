#!/usr/local/bin/Rscript

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(rstan)
library(coda)
library(MCMCglmm)
library(dyndimred)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/pseudogp/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# perform dimreds
dimred_names <- names(dyndimred::list_dimred_methods())
dimred_names <- dimred_names[which(as.logical(params$dimreds))] # 'which()' is to ensure when new dimreds are added, they simply get left out
spaces <- map(dimred_names, ~ dyndimred::dimred(expression, method = ., ndim = 2)) # only 2 dimensions per dimred are allowed

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# fit probabilistic pseudotime model
fit <- pseudogp::fitPseudotime(
  X = spaces,
  smoothing_alpha = params$smoothing_alpha,
  smoothing_beta = params$smoothing_beta,
  iter = params$iter,
  chains = params$chains,
  initialise_from = params$initialise_from,
  pseudotime_var = params$pseudotime_var,
  pseudotime_mean = params$pseudotime_mean
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# extract pseudotime
pst <- rstan::extract(fit, pars = "t")$t
tmcmc <- coda::mcmc(pst)
pseudotime <- MCMCglmm::posterior.mode(tmcmc) %>%
  setNames(rownames(expression))

# return output
output <- lst(
  cell_ids = names(pseudotime),
  pseudotime = pseudotime,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
