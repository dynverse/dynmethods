library(tidyverse)

# generate methods object
methods <- dynwrap::get_ti_methods(ti_packages = "dynmethods") %>%
  select(-method_func)

devtools::use_data(methods, overwrite = TRUE)
devtools::document()
