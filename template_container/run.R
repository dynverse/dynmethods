#!/usr/local/bin/Rscript

task <- dyncli::main()

# load libraries
library(dyncli, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)

library(TEMPLATE)

#####################################
###           LOAD DATA           ###
#####################################
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
