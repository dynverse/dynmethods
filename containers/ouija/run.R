library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(ouija)
library(rstan)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, features_id, iter = 100, response_type = "switch", 
    inference_type = "hmc", normalise_expression = TRUE) 
{
    requireNamespace("ouija")
    requireNamespace("rstan")
    requireNamespace("coda")
    expression <- expression[, features_id]
    rstan::rstan_options(auto_write = TRUE)
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    oui <- ouija::ouija(x = expression, iter = iter, response_type = response_type, 
        inference_type = inference_type, normalise_expression = normalise_expression)
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    pseudotime <- ouija::map_pseudotime(oui) %>% setNames(rownames(expression))
    space <- dimred(expression, method = "pca", ndim = 2)
    k_trace <- rstan::extract(oui$fit, "k")$k
    kmean <- colMeans(k_trace)
    t0 <- rstan::extract(oui$fit, "t0")$t0
    t0_means <- colMeans(t0)
    t0_interval <- coda::HPDinterval(coda::mcmc(t0))
    t0_df <- data_frame(t0_mean = t0_means, lower = t0_interval[, 
        1], upper = t0_interval[, 2], kmean = kmean)
    t0_df$Gene <- colnames(oui$Y[, oui$response_type == "switch"])
    t0_df$Gene <- factor(t0_df$Gene, t0_df$Gene[order(t0_means)])
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_linear_trajectory(pseudotime = pseudotime, t0_df = t0_df) %>% 
        add_dimred(dimred = space) %>% add_timings(timings = tl %>% 
        add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')