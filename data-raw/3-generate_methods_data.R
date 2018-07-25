library(tidyverse)

# generate methods object
methods <- dynwrap::get_ti_methods(as_tibble = FALSE, ti_packages = "dynmethods") %>%
  map(function(method) {
    method$method_func() %>% discard(is.function)
  }) %>%
  dynutils::list_as_tibble()

devtools::use_data(methods, overwrite = TRUE)
devtools::document()
