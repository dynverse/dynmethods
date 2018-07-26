Sys.setenv("R_TESTS" = "")

library(testthat)
library(dynutils)
library(dplyr)
library(purrr)
library(tibble)

test_check("dynmethods")
