library(testthat)
library(dyneval)
library(dynutils)
library(dplyr)
library(ggplot2)
library(purrr)

Sys.setenv("R_TESTS" = "")

test_check("dyneval")

