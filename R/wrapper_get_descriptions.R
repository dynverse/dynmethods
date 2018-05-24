#' Return all TI ti_methods
#'
#' @param as_tibble Whether or not to return the ti_methods as a tibble
#'
#' @importFrom utils lsf.str
#' @export
get_ti_methods <- function(as_tibble = TRUE) {
  requireNamespace("dynmethods")

  ti_methods <- lsf.str(asNamespace("dynmethods"), pattern = "ti_") %>%
    map(~ do.call(., args = list(), envir = asNamespace("dynmethods")))

  if (as_tibble) {
    list_as_tibble(ti_methods)
  } else {
    ti_methods %>% setNames(ti_methods %>% map_chr(~.$short_name))
  }
}
