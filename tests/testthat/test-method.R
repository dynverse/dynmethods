
test_that("A method can be used with docker", {
  dataset <- dynwrap::example_dataset

  config <- babelwhale::create_docker_config()
  babelwhale::set_default_config(config, permanent = FALSE)

  options(dynwrap_backend = "container")

  method <- ti_random()

  trajectory <- infer_trajectory(dataset, method)

  expect_s3_class(trajectory, "dynwrap::with_trajectory")
})