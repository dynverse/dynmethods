context("Testing infer_trajectory")


test_that("Testing infer_trajectory with control methods", {
  data("toy_tasks", package="dyntoy")

  # run with one task and one method
  task <- toy_tasks %>% extract_row_to_list(1)
  method <- ti_comp1()

  model <- infer_trajectory(task, method)
  expect_s3_class(model, "dynwrap::with_trajectory")

  model <- infer_trajectory(task, method, give_priors = c("start_cells"))
  expect_s3_class(model, "dynwrap::with_trajectory")

  expect_error(infer_trajectory(task, method, give_priors = c("to be or not to be")))

  # run with multiple tasks and one method
  models <- infer_trajectories(list(task, task), method)
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)
  expect_setequal(c("task_ix", "method_ix", "model", "method_name", "task_id", "summary"), names(models))

  models <- infer_trajectories(list_as_tibble(list(task, task)), ti_angle())
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)

  # run with multiple methods
  models <- infer_trajectories(task, list(ti_angle(), ti_angle()))
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)

  models <- infer_trajectories(task, list_as_tibble(list(ti_angle(), ti_angle())))
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)

  models <- infer_trajectories(task, c("shuffle", "random"))
  expect_true(is_tibble(models))
  expect_equal(nrow(models), 2)

  expect_message(infer_trajectories(task, c("shoffle")))

  expect_error(infer_trajectories(task, c(1,2,3)))
  expect_error(infer_trajectories(c(1,2,3), c(1,2,3)))

  # run with multiple tasks and multiple methods
  models <- infer_trajectories(
    task = list(task, task, task),
    method = list(ti_angle(), ti_comp1())
  )

  expect_true(is_tibble(models))
  expect_equal(nrow(models), 6)

  # run with multiple tasks and multiple methods with specified parameters
  models <- infer_trajectories(
    task = toy_tasks[c(1,2),],
    method = list_as_tibble(list(ti_angle(), ti_comp1())),
    parameters = list(list(method = "mds"), list(method = "pca"))
  )

  expect_true(is_tibble(models))
  expect_equal(nrow(models), 4)
})
