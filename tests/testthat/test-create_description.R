context("Testing create_description")


test_that("Testing create_description with dummy method", {
  dummy <- dynmethods:::create_description(
    name = "dummy 1",
    short_name = "dum1",
    package_loaded = c("dynverse"),
    package_required = c("tidyverse"),
    par_set = ParamHelpers::makeParamSet(
      ParamHelpers::makeDiscreteParam(id = "param", default = "banana", values = c("apple", "banana", "cherry"))
    ),
    properties = c("space", "trajectory"),
    run_fun = function(counts, param = "fjioiw") param,
    plot_fun = function(out) "cake",
    run_fun_defaults = list(param = "101010")
  )
  expect_equal( dummy$name, "dummy 1" )
  expect_equal( dummy$short_name, "dum1" )
  expect_equal( dummy$package_loaded, "dynverse" )
  expect_equal( dummy$package_required, "tidyverse" )
  expect_is( dummy$par_set, "ParamSet" )
  expect_equal( dummy$properties, c("space", "trajectory") )
  expect_is( dummy$run_fun, "function" )
  # take into account parameter overwriting by parmamset
  expect_equal( dummy$run_fun(NULL), "101010" )
  expect_is( dummy$plot_fun, "function" )
  expect_equal( dummy$plot_fun(NULL), "cake" )
})
