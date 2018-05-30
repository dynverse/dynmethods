context("Testing create_docker_ti_method")

test_that("Testing create_docker_ti_method with compone", {
  method1 <- create_docker_ti_method("dynverse/comp1")
  expect_equal(method1()$short_name, "component_one")

  method2 <- create_docker_ti_method("dynverse/comp1", "test", pull=FALSE)
  expect_equal(method2()$short_name, "test")

  expect_error(create_docker_ti_method("dynverse/comp1", input = "whatever", pull=FALSE))
  expect_error(create_docker_ti_method("dynverse/comp1", output = "whatever", pull=FALSE))
})
