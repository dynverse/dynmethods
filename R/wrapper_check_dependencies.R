#' Check which packages need to be installed for all TI methods to function correctly
#'
#' @export
check_dependencies <- function() {
  for (descr in get_descriptions(as_tibble = FALSE)) {
    required_packages <- c(descr$package_loaded, descr$package_required)
    installed <- required_packages %in% rownames(installed.packages())
    if (any(!installed)) {
      message(sQuote(descr$name), " requires the following packages still to be installed: ", paste(sQuote(required_packages[!installed]), collapse = ", "))
    }
  }
}
