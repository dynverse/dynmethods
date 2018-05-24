#' Get the default parameters of a method
#'
#' @param method A TI method description
#'
#' @export
#'
#' @importFrom testthat expect_true
#' @importFrom ParamHelpers dfRowToList generateDesignOfDefaults
get_default_parameters <- function(method) {
  testthat::expect_true(is_ti_method(method))

  ParamHelpers::dfRowToList(
    ParamHelpers::generateDesignOfDefaults(method$par_set, trafo = TRUE),
    method$par_set,
    i = 1
  )
}
