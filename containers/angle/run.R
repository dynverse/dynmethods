library(dynmethods)
library(dynwrap)
library(jsonlite)
library(readr)

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

method <- do.call(ti_angle, params[names(formals(ti_angle))])
model <- do.call(method$run_fun, data)

write_rds(model, "/output/output.rds")
