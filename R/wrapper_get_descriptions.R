#' Return all TI descriptions
#'
#' @param as_tibble Whether or not to return the descriptions as a tibble
#'
#' @importFrom utils lsf.str
#' @export
get_descriptions <- function(as_tibble = TRUE) {
  requireNamespace("dynmethods")

  descriptions <- lsf.str(asNamespace("dynmethods"), pattern = "description_") %>%
    map(~ do.call(., args = list(), envir = asNamespace("dynmethods")))

  if (as_tibble) {
    list_as_tibble(descriptions)
  } else {
    descriptions %>% setNames(descriptions %>% map_chr(~.$short_name))
  }
}
