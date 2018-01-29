library(dynalysis)
library(tidyverse)

load("~/tmp.RData")

rm(list = setdiff(ls(), c("task", "tasks")))

method <- description_topslam()

def_prms <- ParamHelpers::generateDesignOfDefaults(method$par_set, trafo = TRUE) %>% ParamHelpers::dfRowToList(method$par_set, 1)
list2env(def_prms, environment())


# start_cells <- task$prior_information$start_cells
# counts <- task$counts

out <- execute_method(tasks, method, list())


out[[1]]$summary$error
