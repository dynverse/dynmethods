context("Testing TI method wrappers")

library(ggplot2)

test_that("Descriptions can be retrieved", {
  tib <- get_descriptions()
  expect_that(tib, is_a("tbl"))

  lis <- get_descriptions(as_tibble = FALSE)
  expect_that(lis, is_a("list"))

  for (descr in lis) {
    test_that(paste0("Description ", descr$name), {
      expect_lte(nchar(descr$short_name), 8)
    })
  }
})

methods <- get_descriptions()

for (i in seq_len(nrow(methods))) {
  method <- extract_row_to_list(methods, i)

  test_that(paste0("Checking ", method$short_name), {
    par_set <- method$par_set

    # must be able to generate a 10 random parameters
    design <- ParamHelpers::generateDesign(10, par_set)

    # must be able to generate the default parameters
    design <- ParamHelpers::generateDesignOfDefaults(par_set)

    parset_params <- names(par_set$pars)
    runfun_params <- setdiff(formalArgs(method$run_fun), c("counts", "start_cells", "start_cell", "end_cells", "grouping_assignment", "task"))

    expect_equal( parset_params[parset_params %in% runfun_params], parset_params )
  })
}

test_that("Checking for dependencies does not produce an error", {
  expect_error(check_dependencies(), NA)
})

test_that("Testing create_description with dummy method", {
  dummy <- dyneval:::create_description(
    name = "dummy 1",
    short_name = "dum1",
    package_loaded = c("dynverse"),
    package_required = c("tidyverse"),
    par_set = ParamHelpers::makeParamSet(
      ParamHelpers::makeDiscreteParam(id = "param", default = "banana", values = c("apple", "banana", "cherry"))
    ),
    properties = c("space", "trajectory"),
    run_fun = function(counts, param = "fjioiw") param,
    plot_fun = function(out) "cake",
    override_runfun_params = TRUE
  )
  expect_equal( dummy$name, "dummy 1" )
  expect_equal( dummy$short_name, "dum1" )
  expect_equal( dummy$package_loaded, "dynverse" )
  expect_equal( dummy$package_required, "tidyverse" )
  expect_is( dummy$par_set, "ParamSet" )
  expect_equal( dummy$properties, c("space", "trajectory") )
  expect_is( dummy$run_fun, "function" )
  # take into account parameter overwriting by parmamset
  expect_equal( dummy$run_fun(NULL), "banana" )
  expect_is( dummy$plot_fun, "function" )
  expect_equal( dummy$plot_fun(NULL), "cake" )

  dummy <- dyneval:::create_description(
    name = "dummy 1",
    short_name = "dum1",
    package_loaded = c("dynverse"),
    package_required = c("tidyverse"),
    par_set = ParamHelpers::makeParamSet(
      ParamHelpers::makeDiscreteParam(id = "param", default = "banana", values = c("apple", "banana", "cherry"))
    ),
    properties = c("space", "trajectory"),
    run_fun = function(counts, param = "fjioiw") param,
    plot_fun = function(out) "cake",
    override_runfun_params = FALSE
  )
  # take into account parameter overwriting by parmamset
  expect_equal( dummy$run_fun(NULL), "fjioiw" )
})



test_that("Testing execute_method with dummy method", {
  data("toy_tasks", package="dyntoy")

  dummy <- dyneval:::create_description(
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

      wrap_ti_prediction(
        ti_type = "linear",
        id = "dum2",
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

test_that("Testing timeout of execute_method", {


  timeouter <- dyneval:::create_description(
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

      wrap_ti_prediction(
        ti_type = "linear",
        id = "dum3",
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
