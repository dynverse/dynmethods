#' Return all TI descriptions
#'
#' @param as_tibble Whether or not to return the descriptions as a tibble
#'
#' @importFrom utils lsf.str
#' @export
get_descriptions <- function(as_tibble = TRUE) {
  requireNamespace("dynmethods")
  functions <- lsf.str(asNamespace("dynmethods"))
  description_functions <- functions[grep("description_", functions)]
  descriptions <- lapply(description_functions, function(fun_name) {
    do.call(fun_name, args = list(), envir = asNamespace("dynmethods"))
  })
  if (as_tibble) {
    list_as_tibble(descriptions)
  } else {
    descriptions %>% setNames(descriptions %>% map_chr(~.$short_name))
  }
}
