context("Testing execute_method")


test_that("Testing execute_method with dummy method", {
  data("toy_tasks", package="dyntoy")

  dummy <- dynmethods:::create_description(
    name = "dummy 2",
    short_name = "dum2",
    package_loaded = c("dplyr"),
    package_required = c(),
    par_set = ParamHelpers::makeParamSet(
      ParamHelpers::makeDiscreteParam(id = "aggr_fun", values = c("mean", "median"), default = "mean")
    ),
    properties = c(),
    run_fun = function(counts, aggr_fun = "mean") {
      pt <- apply(counts, 1, aggr_fun) %>% dynutils::scale_minmax()

      milestone_ids <- c("start", "end")
      milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1, directed=TRUE)
      progressions <- data_frame(cell_id = names(pt), from = milestone_ids[[1]], to = milestone_ids[[2]], percentage = pt)

      wrap_prediction_model(
        cell_ids = rownames(counts),
        milestone_ids = milestone_ids,
        milestone_network = milestone_network,
        progressions = progressions
      )
    },
    plot_fun = function(out) {
      ggplot(out$progressions, aes(percentage, percentage)) +
        geom_point(size = 2.2) +
        geom_point(aes(colour = percentage), size = 2) +
        scale_colour_distiller(palette = "RdBu")
    }
  )

  method_outs <- execute_method(toy_tasks, dummy, parameters = list(aggr_fun = "median"), timeout = 1e6)

  for (i in seq_along(method_outs)) {
    method_out <- method_outs[[i]]

    expect_true( dynutils::is_ti_data_wrapper(method_out$model) )
    expect_is( method_out$summary, "data.frame" )

    pdf("/dev/null")
    expect_error( print(dummy$plot_fun(method_out$model)), NA)
    dev.off()

    expect_equal( nrow(method_out$summary), 1 )
  }
})


test_that("Testing prior passing for execute_method", {
  data("toy_tasks", package="dyntoy")

  dummy <- dynmethods:::create_description(
    name = "dummy 2",
    short_name = "dum2",
    package_loaded = c("dplyr"),
    package_required = c(),
    par_set = ParamHelpers::makeParamSet(
      ParamHelpers::makeDiscreteParam(id = "aggr_fun", values = c("mean", "median"), default = "mean")
    ),
    properties = c(),
    run_fun = function(
      expression,
      start_milestones,
      start_cells,
      end_milestones,
      end_cells,
      grouping_assignment,
      grouping_network,
      marker_feature_ids,
      n_branches,
      time,
      timecourse,
      aggr_fun = "mean"
    ) {
      pt <- apply(expression, 1, aggr_fun) %>% dynutils::scale_minmax()

      milestone_ids <- c("start", "end")
      milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1, directed=TRUE)
      progressions <- data_frame(cell_id = names(pt), from = milestone_ids[[1]], to = milestone_ids[[2]], percentage = pt)

      wrap_prediction_model(
        milestone_ids = milestone_ids,
        milestone_network = milestone_network,
        progressions = progressions,
        cell_ids = names(pt)
      )
    },
    plot_fun = function(out) {
      ggplot(out$progressions, aes(percentage, percentage)) +
        geom_point(size = 2.2) +
        geom_point(aes(colour = percentage), size = 2) +
        scale_colour_distiller(palette = "RdBu")
    }
  )

  method_outs <- execute_method(toy_tasks, dummy, parameters = list(aggr_fun = "median"), timeout = 1e6)

  for (i in seq_along(method_outs)) {
    method_out <- method_outs[[i]]

    expect_true( dynutils::is_ti_data_wrapper(method_out$model) )
    expect_is( method_out$summary, "data.frame" )

    pdf("/dev/null")
    expect_error( print(dummy$plot_fun(method_out$model)), NA)
    dev.off()

    expect_equal( nrow(method_out$summary), 1 )
  }
})

test_that("Testing timeout functionality of execute_method with dummy wrapper", {
  timeouter <- dynmethods:::create_description(
    name = "timeouter",
    short_name = "timeout",
    package_loaded = c("dplyr"),
    package_required = c(),
    par_set = ParamHelpers::makeParamSet(
      ParamHelpers::makeNumericParam(id = "sleep_time", lower = 1, upper = 20, default = 10)
    ),
    properties = c(),
    run_fun = function(counts, sleep_time = 10) {
      Sys.sleep(sleep_time)

      pt <- apply(counts, 1, "mean") %>% dynutils::scale_minmax()
      milestone_ids <- c("start", "end")
      milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1, directed=TRUE)
      progressions <- data_frame(cell_id = names(pt), from = milestone_ids[[1]], to = milestone_ids[[2]], percentage = pt)

      wrap_prediction_model(
        cell_ids = rownames(counts),
        milestone_network = milestone_network,
        progressions = progressions,
        milestone_ids = milestone_ids
      )
    },
    plot_fun = function(out) {
      ggplot(out$progressions, aes(percentage, percentage)) +
        geom_point(size = 2.2) +
        geom_point(aes(colour = percentage), size = 2) +
        scale_colour_distiller(palette = "RdBu")
    }
  )

  data("toy_tasks", package="dyntoy")
  toy_tasks <- toy_tasks[1:6,]

  num_datasets <- nrow(toy_tasks)

  # should take about 20 seconds, but will timeout after 10
  out <- execute_method(toy_tasks, timeouter, parameters = list(sleep_time = 10), timeout = 1)
  models <- out %>% map_df(~ .$model)
  summaries <- out %>% map_df(~.$summary)
  expect_equal(nrow(models), 0)
  expect_equal(nrow(summaries), num_datasets)
  for (i in seq_len(nrow(summaries))) {
    expect_false(is.null(summaries$error[[i]]))
    expect_is(summaries$error[[i]], "error")
  }

  # this should finish correctly
  method_outs <- execute_method(toy_tasks, timeouter, parameters = list(sleep_time = 1), timeout = 10)

  for (i in seq_along(method_outs)) {
    method_out <- method_outs[[i]]

    expect_true( dynutils::is_ti_data_wrapper(method_out$model) )
    expect_is( method_out$summary, "data.frame" )

    pdf("/dev/null")
    expect_error( print(timeouter$plot_fun(method_out$model)), NA)
    dev.off()

    expect_equal( nrow(method_out$summary), 1 )
    expect_true( is.null(method_out$summary$error[[1]]) )
  }
})

tasks <- dyntoy::generate_toy_datasets(trajectory_types = c("simple_linear"), num_replicates = 1, num_cells = 100, num_genes = 51)
# tasks <- dyntoy::toy_tasks %>% filter(trajectory_type == "linear") %>% slice(1)
methods <- get_descriptions(as_tibble = FALSE)
for (method in methods) {
  test_that(pritt("Testing whether {method$short_name} is able to run on simple data"), {
    params <- ParamHelpers::generateDesignOfDefaults(method$par_set, trafo = TRUE) %>% ParamHelpers::dfRowToList(method$par_set, 1)
    out <- execute_method(tasks, method, parameters = params, timeout = 100)
    error <- out[[1]]$summary$error[[1]]
    error
    expect_null(error)
    # if erroring:
    # list2env(params, globalenv())
    # counts <- toy$counts[[1]]
    # expression <- toy$expression[[1]]
  })
}

# TODO: implement test for checking whether zerod columns work
