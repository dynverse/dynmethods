context("Testing check_dependencies")

test_that("Checking for dependencies does not produce an error", {
  expect_error(check_dependencies(), NA)
})
