context("Testing get_ti_methods")

# only run this test on the maintainer's development environment
skip_if_not(Sys.info()[["user"]] %in% c("rcannood", "wouters"))

methods <- get_ti_methods(ti_packages = "dynmethods")

for (i in seq_len(nrow(methods))) {
  method <- extract_row_to_list(methods, i)$fun()

  test_that(pritt("Checking whether {method$short_name} can generate parameters"), {
    par_set <- method$par_set

    # must be able to generate a 3 random parameters
    design <- ParamHelpers::generateDesign(3, par_set)

    # must be able to generate the default parameters
    design <- ParamHelpers::generateDesignOfDefaults(par_set)

    parset_params <- names(par_set$pars)
    runfun_params <- setdiff(formalArgs(method$run_fun), c("counts", "start_id", "start_cell", "end_id", "groups_id", "dataset"))

    expect_equal( parset_params[parset_params %in% runfun_params], parset_params )
  })
}
