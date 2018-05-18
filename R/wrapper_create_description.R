#' Create a TI method description
#'
#' @param name The name of the TI method
#' @param short_name A short name for the method, max 8 characters
#' @param package_loaded The packages that need to be loaded before executing the method
#' @param package_required The packages that need to be installed before executing the method
#' @param par_set A bunch of parameters created by [ParamHelpers::makeParamSet()]
#' @param properties Several descriptive properties of the method. WIP.
#' @param run_fun A function to run the TI, needs to have 'counts' as its first param.
#' @param plot_fun A function to plot the results of a TI, needs to have 'prediction' as its first param.
#' @param override_runfun_params Whether or not to override the default parameters
#' of `run_fun` with those described in `par_set`.
create_description <- function(
  name,
  short_name,
  package_loaded,
  package_required,
  par_set,
  properties,
  run_fun,
  plot_fun,
  override_runfun_params = TRUE
) {
  if (override_runfun_params) {
    default_params <- par_set %>%
      generateDesignOfDefaults(trafo = TRUE) %>%
      ParamHelpers::dfRowToList(par_set, 1)

    if(!all(names(default_params) %in% formalArgs(run_fun))) {
      stop("Not all default params described in par_set are listed in the run_fun.")
    }

    formals(run_fun)[names(default_params)] <- default_params
  }
  desc <- lst(
    name,
    short_name,
    package_loaded,
    package_required,
    par_set,
    properties,
    run_fun,
    plot_fun
  )
  class(desc) <- c("dynmethod::description", class(desc))
  desc
}

#' Tests whether an object is a TI method description
#'
#' @param object The object to be tested
#'
#' @export
is_description <- function(object) {
  "dynmethod::description" %in% class(object)
}
