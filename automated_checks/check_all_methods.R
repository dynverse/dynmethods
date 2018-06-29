library(dynmethods)
library(tidyverse)
data(methods, package = "dynmethods")

methods <- methods %>% filter(containerised)

walk(
  methods$method_id, function(method_id) {
    model <- NULL
    file.remove(list.files("automated_checks/", pattern = paste0(method_id, "\\..*"), full.names = TRUE))

    rmarkdown::render(
      "automated_checks/check_method.Rmd",
      output_format = rmarkdown::github_document(),
      params = list(method_id = method_id),
      knit_root_dir = "./",
      output_file = paste0("", method_id, ".md"),
      quiet = TRUE
    )

    rm(list=ls())
  }
)

