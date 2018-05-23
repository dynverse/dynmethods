context("Testing all wrappers")


# Todo: re-enable this
# tasks <- dyntoy::toy_tasks %>% filter(model == "simple_linear") %>% slice(1)
# methods <- get_descriptions(as_tibble = FALSE)
# for (method in methods) {
#   test_that(pritt("Testing whether {method$short_name} is able to run on simple data"), {
#     params <- ParamHelpers::generateDesignOfDefaults(method$par_set, trafo = TRUE) %>% ParamHelpers::dfRowToList(method$par_set, 1)
#     out <- infer_trajectory(tasks, method, parameters = params)
#     error <- out[[1]]$summary$error[[1]]
#     error
#     expect_null(error)
#     # if erroring:
#     # list2env(params, globalenv())
#     # counts <- toy$counts[[1]]
#     # expression <- toy$expression[[1]]
#   })
# }
