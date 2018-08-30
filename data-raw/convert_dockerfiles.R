library(tidyverse)
library(furrr)
library(dynmethods)
library(dynwrap)

files <- list.files("containers", pattern = "Dockerfile", recursive = TRUE, full.names = TRUE)

# create singularity scripts from dockerfiles
walk(files, function(file) {
  new_file <- gsub("containers/([^/]*)/.*", "containers/\\1/Singularity.\\1", file)
  dynwrap:::.container_dockerfile_to_singularityrecipe(file, new_file)
})
