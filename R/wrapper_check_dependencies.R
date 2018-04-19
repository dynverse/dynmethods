#' Check which packages need to be installed for all TI methods to function correctly
#'
#' @export
check_dependencies <- function() {
  for (descr in get_descriptions(as_tibble = FALSE)) {
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
  installed <- required_packages %in% rownames(installed.packages())
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

  remotes <- desc::desc_get_remotes(find.package("dynmethods")) %>%
    set_names(., stringr::str_replace(., ".*/([:alpha:]*).*", "\\1"))

  for (dependency in dependencies[dependencies %in% names(remotes)]) {
    devtools::install_github(remotes[[dependency]])
  }
  for (dependency in dependencies[!dependencies %in% names(remotes)]) {
    devtools::install_cran(dependency)
  }

  message("Installed ", paste0(dependencies, collapse = ", "))
}

