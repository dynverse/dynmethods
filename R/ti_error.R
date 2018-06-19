#' Inferring trajectories with Control: error
#'
#' This control method will always produce an error.
#'
#' @param dummy_param This parameter does not do anything.
#'
#' @export
ti_error <- create_ti_method(
  name = "Control: error",
  short_name = "error",
  package_loaded = c(),
  package_required = c(),
  trajectory_types = c("linear", "bifurcation", "convergence"),
  type = "control",
  authors = list(
    list(
      given = "Robrecht",
      family = "Cannoodt",
      email = "rcannood@gmail.com",
      ORCID = "0000-0003-3641-729X",
      github = "rcannood"
    ),
    list(
      given = "Wouter",
      family = "Saelens",
      email = "wouter.saelens@ugent.be",
      ORCID = "0000-0002-7114-6248",
      github = "zouter"
    )
  ),
  par_set = makeParamSet(
    makeNumericParam(id = "dummy_param", lower = 0, default = 0.5, upper = 1)
  ),
  run_fun = "dynmethods::run_error",
  plot_fun = dynplot::plot_default
)

run_error <- function(counts, dummy_param) {
  stop("This control method always errors.")
}

