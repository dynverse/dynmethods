library(dynwrap)
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

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, dimreds = c(TRUE, TRUE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE), chains = 3, iter = 100, smoothing_alpha = 10, 
    smoothing_beta = 3, pseudotime_mean = 0.5, pseudotime_var = 1, 
    initialise_from = "random") 
{
    requireNamespace("pseudogp")
    requireNamespace("rstan")
    requireNamespace("coda")
    requireNamespace("MCMCglmm")
    dimred_names <- names(dyndimred::list_dimred_methods())[as.logical(dimreds)]
    spaces <- map(dimred_names, ~dimred(expression, method = ., 
        ndim = 2))
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    fit <- pseudogp::fitPseudotime(X = spaces, smoothing_alpha = smoothing_alpha, 
        smoothing_beta = smoothing_beta, iter = iter, chains = chains, 
        initialise_from = initialise_from, pseudotime_var = pseudotime_var, 
        pseudotime_mean = pseudotime_mean)
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    pst <- rstan::extract(fit, pars = "t")$t
    tmcmc <- coda::mcmc(pst)
    pseudotime <- MCMCglmm::posterior.mode(tmcmc) %>% setNames(rownames(expression))
    pst <- rstan::extract(fit, pars = "t", permute = FALSE)
    lambda <- rstan::extract(fit, pars = "lambda", permute = FALSE)
    sigma <- rstan::extract(fit, pars = "sigma", permute = FALSE)
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_linear_trajectory(pseudotime = pseudotime, spaces = spaces, 
            chains = chains, pst = pst, lambda = lambda, sigma = sigma) %>% 
        add_timings(timings = tl %>% add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')