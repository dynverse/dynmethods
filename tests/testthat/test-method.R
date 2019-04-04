context("Test one methods")

test_that("A method can be used with docker", {
  dataset <- dynwrap::example_dataset

  babelwhale::test_docker_installation(detailed = TRUE)

  config <- babelwhale::create_docker_config()
  babelwhale::set_default_config(config)

  method <- ti_comp1()

  trajectory <- infer_trajectory(dataset, method)
})

# test_that("A method can be used with singularity", {
#   dataset <- dynwrap::example_dataset
#
#   babelwhale::test_singularity_installation(detailed = TRUE)
#
#   config <- babelwhale::create_singularity_config()
#   babelwhale::set_default_config(config)
#
#   method <- ti_comp1()
#
#   trajectory <- infer_trajectory(dataset, method)
# })
