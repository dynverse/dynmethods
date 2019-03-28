context("Test whether methods-object exists")

test_that("methods-object exists", {
  expect_is(dynmethods::methods$method_id, "character")
})
