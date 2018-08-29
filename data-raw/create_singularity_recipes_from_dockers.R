library(tidyverse)
library(furrr)
library(dynmethods)
library(dynwrap)

files <- list.files("containers", pattern = "Dockerfile", recursive = TRUE, full.names = TRUE)

# create singularity scripts from dockerfiles
walk(files, function(file) {
  new_file <- gsub("containers/([^/]*)/.*", "containers/\\1/Singularity.\\1", file)
  lines <-
    readr::read_lines(file) %>%
    discard(~ grepl("^ *$", .))

  from <- str_subset(lines, "^FROM ") %>% str_replace_all("FROM ", "From: ")
  environment <- str_subset(lines, "^ENV ") %>% str_replace_all("ENV ", "    export ")
  labels <- str_subset(lines, "^LABEL ") %>% str_replace_all("LABEL ", "    ")
  post <- str_subset(lines, "^RUN ") %>% str_replace_all("RUN ", "    ")
  files <- str_subset(lines, "^ADD ") %>% str_replace_all("ADD ", "    ")
  runscript <- str_subset(lines, "^ENTRYPOINT ") %>% str_replace_all("ENTRYPOINT ", "    exec ")

  sing_lines <- c(
    "Bootstrap: shub",
    "",
    from,
    "",
    "%environment",
    environment,
    "",
    "%labels",
    labels,
    "",
    "%post",
    post,
    "",
    "%files",
    files,
    "",
    "%runscript",
    runscript
  )

  readr::write_lines(sing_lines, new_file)
})
