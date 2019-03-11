#!/usr/local/bin/Rscript

library(dyncli)

#####################################
###           LOAD DATA           ###
#####################################

# load data
task <- dyncli::main()

# load libraries
library(dynwrap)
library(dplyr)
library(purrr)

library(TEMPLATE)

expression <- task$expression
params <- task$params
priors <- task$priors

# TIMING: done with preproc
timings <- list(method_afterpreproc = Sys.time())

#####################################
###        INFER TRAJECTORY       ###
#####################################

# do TI calculations

# TIMING: done with trajectory inference
timings$method_aftermethod <- Sys.time()

#####################################
###     SAVE OUTPUT TRAJECTORY    ###
#####################################


# save output
output <-
  wrap_data(
    cell_ids = rownames(expression)
  ) %>%
  add_dimred(
    dimred = dimred
  )%>%
  add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep
  )  %>%
  add_timings(
    timings = timings
  )

dyncli::write_output(output, task$output)
