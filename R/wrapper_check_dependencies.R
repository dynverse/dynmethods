#' Check which packages need to be installed for all TI methods to function correctly
#'
#' @export
check_dependencies <- function() {
  for (descr in get_ti_methods(as_tibble = FALSE)) {
    check_dependency(descr)
  }
}

#' Check which packages need to be installed for a TI method to function correctly
#'
#' @param descr Description of the method
#'
#' @export
check_dependency <- function(descr) {
  required_packages <- c(descr$package_loaded, descr$package_required)
  installed <- check_packages(required_packages)

  if (any(!installed)) {
    message(sQuote(descr$name), " requires the following packages still to be installed: ", paste(sQuote(required_packages[!installed]), collapse = ", "))
  }

  required_packages
}


#' Install dependencies for a method
#'
#' @inheritParams check_dependency
#'
#' @export
install_dependencies <- function(descr) {
  dependencies <- check_dependency(descr)

  install_packages(dependencies, "dynmethods")
}
