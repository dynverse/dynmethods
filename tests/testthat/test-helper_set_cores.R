context("Testing helper_set_cores")

test_that("Testing helper_set_cores", {
  dynmethods:::set_cores("111")
  expect_equal(Sys.getenv("MKL_NUM_THREADS"), "111")
  expect_equal(Sys.getenv("NUMEXPR_NUM_THREADS"), "111")
  expect_equal(Sys.getenv("OMP_NUM_THREADS"), "111")

  dynmethods:::set_cores("123")
  expect_equal(Sys.getenv("MKL_NUM_THREADS"), "123")
  expect_equal(Sys.getenv("NUMEXPR_NUM_THREADS"), "123")
  expect_equal(Sys.getenv("OMP_NUM_THREADS"), "123")
})
