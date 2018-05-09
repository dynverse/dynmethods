library(testthat)
library(dynmethods)
library(dynutils)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)

Sys.setenv("R_TESTS" = "")

test_check("dynmethods")

