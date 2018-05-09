context("Testing infer_trajectory")


test_that("Testing infer_trajectory with control methods", {
  data("toy_tasks", package="dyntoy")

  task <- toy_tasks %>% extract_row_to_list(1)

  # run with one task and one method
  model <- infer_trajectory(task, description_angle())
  expect_s3_class(model, "dynwrap::with_trajectory")

  # run with multiple tasks and one method
  models <- infer_trajectories(list(task, task), description_angle())
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)

  models <- infer_trajectories(list_as_tibble(list(task, task)), description_angle())
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)

  # run with multiple methods
  models <- infer_trajectories(task, list(description_angle(), description_angle()))
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)

  models <- infer_trajectories(task, list_as_tibble(list(description_angle(), description_angle())))
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)

  models <- infer_trajectories(task, c("angle", "angle"))
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)
})
