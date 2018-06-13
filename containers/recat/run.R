library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(reCAT)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, num_cores = 1, TSPFold = 2, beginNum = 10, 
    endNum = 15, step_size = 2, base_cycle_range_start = 6, base_cycle_range_end = 9, 
    max_num = 300, clustMethod = "GMM") 
{
    requireNamespace("reCAT")
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    result <- reCAT::bestEnsembleComplexTSP(test_exp = expression, 
        TSPFold = TSPFold, beginNum = beginNum, endNum = endNum, 
        base_cycle_range = base_cycle_range_start:base_cycle_range_end, 
        step_size = step_size, max_num = max_num, clustMethod = clustMethod, 
        threads = num_cores, output = FALSE)
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    pseudotime <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], 
        ] %>% set_names(rownames(expression))
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_cyclic_trajectory(pseudotime = pseudotime) %>% add_timings(timings = tl %>% 
        add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')