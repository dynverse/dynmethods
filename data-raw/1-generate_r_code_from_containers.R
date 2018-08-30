library(tidyverse)
library(dynwrap)
library(furrr)

source("data-raw/1a-helper_functions.R")

files <- list.files("containers", pattern = "definition.yml", recursive = TRUE, full.names = TRUE)

# iterate over the containers and generate R scripts for each of them
definitions <-
  map(files, function(file) {
    cat(file, "\n", sep = "")

    # read the docker repository
    repo <- yaml::read_yaml(file)$docker_repository

    # fetch definition /with/ digests
    definition <- dynwrap:::.container_get_definition(repo)

    # generate file from definition
    generate_file_from_container(definition)

    # return the definition
    definition
  })
methods <- dynutils::list_as_tibble(definitions)
method_versions <- methods %>%
  select(docker_repository, version) %>%
  deframe()

devtools::use_data(methods, overwrite = TRUE)
devtools::use_data(method_versions, overwrite = TRUE)

# don't forget to regenerate the documentation
devtools::document()
devtools::install(dependencies = FALSE)
