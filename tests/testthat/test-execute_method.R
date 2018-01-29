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

      # TIMING: done with preproc
      tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

      milestone_ids <- c("start", "end")
      milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1, directed=TRUE)
      progressions <- data_frame(cell_id = names(pt), from = milestone_ids[[1]], to = milestone_ids[[2]], percentage = pt)

      # TIMING: done with method
      tl <- tl %>% add_timing_checkpoint("method_aftermethod")

      # TIMING: after postproc
      tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

      wrap_prediction_model(
        cell_ids = rownames(counts),
        milestone_ids = milestone_ids,
        milestone_network = milestone_network,
        progressions = progressions
      ) %>% attach_timings_attribute(tl)
    },
    plot_fun = function(out) {
      ggplot(out$progressions, aes(percentage, percentage)) +
        geom_point(size = 2.2) +
        geom_point(aes(colour = percentage), size = 2) +
        scale_colour_distiller(palette = "RdBu")
    }
  )

  method_outs <- execute_method(
    tasks = toy_tasks,
    method = dummy,
    parameters = list(aggr_fun = "median")
  )

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

      # TIMING: done with preproc
      tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

      milestone_ids <- c("start", "end")
      milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1, directed=TRUE)
      progressions <- data_frame(cell_id = names(pt), from = milestone_ids[[1]], to = milestone_ids[[2]], percentage = pt)

      # TIMING: done with method
      tl <- tl %>% add_timing_checkpoint("method_aftermethod")

      # TIMING: after postproc
      tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

      wrap_prediction_model(
        milestone_ids = milestone_ids,
        milestone_network = milestone_network,
        progressions = progressions,
        cell_ids = names(pt)
      ) %>% attach_timings_attribute(tl)
    },
    plot_fun = function(out) {
      ggplot(out$progressions, aes(percentage, percentage)) +
        geom_point(size = 2.2) +
        geom_point(aes(colour = percentage), size = 2) +
        scale_colour_distiller(palette = "RdBu")
    }
  )

  method_outs <- execute_method(toy_tasks, dummy, parameters = list(aggr_fun = "median"))

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


tasks <- dyntoy::toy_tasks %>% filter(model == "simple_linear") %>% slice(1)
methods <- get_descriptions(as_tibble = FALSE)
for (method in methods) {
  test_that(pritt("Testing whether {method$short_name} is able to run on simple data"), {
    params <- ParamHelpers::generateDesignOfDefaults(method$par_set, trafo = TRUE) %>% ParamHelpers::dfRowToList(method$par_set, 1)
    out <- execute_method(tasks, method, parameters = params)
    error <- out[[1]]$summary$error[[1]]
    error
    expect_null(error)
    # if erroring:
    # list2env(params, globalenv())
    # counts <- toy$counts[[1]]
    # expression <- toy$expression[[1]]
  })
}
