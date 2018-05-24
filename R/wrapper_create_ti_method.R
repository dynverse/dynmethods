#' Create a TI method description
#'
#' @param name The name of the TI method
#' @param short_name A short name for the method, max 8 characters
#' @param package_loaded The packages that need to be loaded before executing the method
#' @param package_required The packages that need to be installed before executing the method
#' @param par_set A bunch of parameters created by [ParamHelpers::makeParamSet()]
#' @param run_fun A function to run the TI, needs to have 'counts' as its first param.
#' @param plot_fun A function to plot the results of a TI, needs to have 'prediction' as its first param.
#'   of `run_fun` with those described in `par_set`.
create_ti_method <- function(
  name,
  short_name,
  package_loaded,
  package_required,
  par_set,
  run_fun,
  plot_fun,
  properties = NA
) {

  desc <- lst(
    name,
    short_name,
    package_loaded,
    package_required,
    par_set
  ) %>% add_class("dynmethod::ti_method")

  default_params <- par_set %>%
    generateDesignOfDefaults(trafo = TRUE) %>%
    ParamHelpers::dfRowToList(par_set, 1)

  ti_fun_constructor_with_params <- function(...) {
    run_fun <- get_function(run_fun)

    plot_fun <- get_function(plot_fun)

    # get the parameters from this function
    run_fun_defaults <- as.list(environment())[formalArgs(ti_fun_constructor_with_params)]

    # override default parameters in the run_fun
    formals(run_fun)[names(run_fun_defaults)] <- run_fun_defaults

    # supply run_fun to the description
    desc$run_fun <- run_fun
    desc$plot_fun <- plot_fun

    # return the description
    desc
  }

  formals(ti_fun_constructor_with_params) <- default_params

  ti_fun_constructor_with_params
}

#' Tests whether an object is a TI method description
#'
#' @param object The object to be tested
#'
#' @export
is_ti_method <- function(object) {
  "dynmethod::ti_method" %in% class(object)
}

get_function <- function(fun) {
  if (is.character(fun)) {
    if (grepl("::", fun)) {
      parts <- strsplit(fun, "::")[[1]]
      get(parts[[2]], envir = asNamespace(parts[[1]]))
    } else {
      get(fun)
    }
  } else if (is.function(fun)) {
    fun
  } else {
    stop(sQuote("fun"), " has an incompatible format.")
  }
}
