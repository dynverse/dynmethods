library(tidyverse)
library(furrr)
library(dynmethods)
library(dynwrap)

bump_version_numbers <- function (by) {
  files <- list.files("../methods", pattern = "Dockerfile", recursive = TRUE, full.names = TRUE)

  increment <- function(v) {
    function(x) {
      xv <- x %>%
        str_replace_all("[^0-9\\.]", "") %>%
        strsplit(".", fixed = TRUE) %>%
        first() %>%
        as.integer()
      paste0("version ", paste0(xv + v, collapse = "."))
    }
  }

  walk(files, function(file) {
    readr::read_lines(file) %>%
      str_replace_all("version .*", increment(c(0, 0, 1))) %>%
      readr::write_lines(file)

    new_file <- gsub("Dockerfile", "Singularity", file)
    dynwrap:::.container_dockerfile_to_singularityrecipe(file, new_file)
  })
}

bump_version_numbers(c(0, 0, 1))
