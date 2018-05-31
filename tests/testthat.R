Sys.setenv("R_TESTS" = "")

library(testthat)
library(dynmethods)
library(dynutils)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)

test_check("dynmethods")
