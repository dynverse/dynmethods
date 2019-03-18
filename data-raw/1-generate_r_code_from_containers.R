library(tidyverse)
library(dynwrap)
library(dynutils)
library(furrr)

source("data-raw/1a-helper_functions.R")

files <- list.files("../methods/", pattern = "definition.yml", recursive = TRUE, full.names = TRUE)

# iterate over the containers and generate R scripts for each of them
definitions <-
  map(files, function(file) {
    cat(file, "\n", sep = "")

    definition <- create_ti_method_definition(definition = file, script = NULL, return_function = FALSE)
    version <- file %>% str_replace("definition.yml", "version") %>% read_lines() %>% str_replace("VERSION=", "")

    # generate file from definition
    generate_file_from_container(definition, version)

    definition$version <- version

    # return the definition
    definition
  })
methods <- dynutils::list_as_tibble(definitions)

for (n in rev(c("method", "wrapper", "container", "manuscript"))) {
  cat("Processing ", n, "\n", sep = "")
  newtib <-
    methods[[n]] %>%
    map(function(x) { if (is.list(x)) x else list() }) %>%
    list_as_tibble()
  newtib[[".object_class"]] <- NULL

  colnames(newtib) <- paste0(n, "_", colnames(newtib))

  methods <- bind_cols(newtib, methods)
  methods[[n]] <- NULL
}
methods[c(".object_class", "run")] <- NULL

usethis::use_data(methods, overwrite = TRUE)

# don't forget to regenerate the documentation
devtools::document()
devtools::install(dependencies = FALSE)
