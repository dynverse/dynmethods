context("Testing create_ti_method")


test_that("Testing create_ti_method with dummy method", {
  dummy <- dynmethods:::create_ti_method(
    name = "dummy 1",
    short_name = "dum1",
    package_loaded = c("dynverse"),
    package_required = c("tidyverse"),
    par_set = ParamHelpers::makeParamSet(
      ParamHelpers::makeDiscreteParam(id = "param", default = "banana", values = c("apple", "banana", "cherry"))
    ),
    run_fun = function(counts, param = "fjioiw") param,
    plot_fun = function(out) "cake"
  )

  dummy_instance <- dummy()

  expect_equal( dummy_instance$name, "dummy 1" )
  expect_equal( dummy_instance$short_name, "dum1" )
  expect_equal( dummy_instance$package_loaded, "dynverse" )
  expect_equal( dummy_instance$package_required, "tidyverse" )
  expect_is( dummy_instance$par_set, "ParamSet" )
  expect_is( dummy_instance$run_fun, "function" )
  # take into account parameter overwriting by parmamset
  expect_equal( dummy_instance$run_fun(NULL), "banana" )
  expect_is( dummy_instance$plot_fun, "function" )
  expect_equal( dummy_instance$plot_fun(NULL), "cake" )

  dummy_instance2 <- dummy(param = "101010")
  expect_equal( dummy_instance2$run_fun(NULL), "101010" )
})
