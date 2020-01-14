library(tidyverse)
library(dynwrap)
library(dynutils)
library(furrr)

source("data-raw/1a-helper_functions.R")

files <- list.files("../methods/", pattern = "Dockerfile", recursive = TRUE, full.names = TRUE)

# iterate over the containers and generate R scripts for each of them
# this loads in the current version from the version files
definitions <-
  map(files, function(docker_file) {
    root_folder <- paste0(dirname(docker_file), "/")
    package_dir <- paste0(root_folder, "package/")
    definition_file <- paste0(root_folder, "definition.yml")
    remotes_file <- paste0(root_folder, "remotes.yml")
    version_file <- paste0(root_folder, "version")

    cat("Processing ", root_folder, "\n", sep = "")

    definition <-
      if (file.exists(definition_file)) {
        create_ti_method_definition(definition = definition_file, script = NULL, return_function = FALSE)

      } else if (file.exists(remotes_file)) {
        remotes <- yaml::read_yaml(remotes_file)

        getFromNamespace(remotes$function_name, remotes$name)()
      }

    # use package version if specified
    if (file.exists(package_dir)) {
      version <- desc::desc(package_dir)$get_version()
    } else {
      version <- version_file %>% read_lines() %>% str_replace("VERSION=", "")
    }

    # generate file from definition
    generate_file_from_container(definition, version)

    definition$version <- version

    # return the definition
    definition
  })
methods <- dynutils::list_as_tibble(definitions)

for (n in rev(c("method", "wrapper", "container", "manuscript", "package"))) {
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
