library(tidyverse)

# generate methods object
methods <- dynwrap::get_ti_methods(as_tibble = FALSE) %>%
  map(function(method) {
    method$method_func() %>% discard(is.function)
  }) %>%
  dynutils::list_as_tibble()

usethis::use_data(methods, overwrite = TRUE)
