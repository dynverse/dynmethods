context("Testing wrap_prediction_model")

cell_ids <- letters[1:10]
cell_info <- data_frame(cell_id = cell_ids, banana = 2, mango = runif(length(cell_ids)))
peach <- 1
tomato <- list(seeds = 102, juice = "over 9000!")

test_that("Testing wrap_prediction_model", {
  out <- wrap_prediction_model(
    cell_ids = cell_ids,
    cell_info = cell_info,
    peach = peach,
    tomato = tomato
  )

  expect_equal(out$cell_ids, cell_ids)
  expect_equal(out$cell_info, cell_info)
  expect_equal(out$peach, peach)
  expect_equal(out$tomato, tomato)

  expect_true(is_prediction(out))
  expect_false(is_prediction(list(1,2,4,5)))
})
