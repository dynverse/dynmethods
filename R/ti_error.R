#' Description for error
#'
#' @importFrom dynplot plot_default
#'
#' @export
description_error <- function() create_description(
  name = "Control: error",
  short_name = "error",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeNumericParam(id = "dummy_param", lower = 0, default = 0.5, upper = 1)
  ),
  properties = c(),
  run_fun = run_error,
  plot_fun = dynplot::plot_default
)

run_error <- function(counts, dummy_param = .5) {
  stop("This control method always errors.")
}

