library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(Waterfall)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, num_clusters = 10) 
{
    requireNamespace("Waterfall")
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    ps <- Waterfall::pseudotimeprog.foo(t(expression), k = num_clusters)
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_linear_trajectory(pseudotime = ps$pseudotime %>% 
            setNames(rownames(expression)), ps = ps) %>% add_timings(timings = tl %>% 
        add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')