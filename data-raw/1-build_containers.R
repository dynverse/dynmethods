library(tidyverse)
library(furrr)

plan(multiprocess)

definition_files <- list.files("containers", pattern = "definition.yml", recursive = TRUE, full.names = TRUE)

future_map(definition_files, function(file) {
  definition <- yaml::read_yaml(file)
  docker_repo <- definition$docker_repository
  folder <- str_replace(file, "/definition.yml$", "")

  system(paste0("docker pull ", docker_repo))
  system(paste0("docker build ", folder, " -t ", docker_repo))
  system(paste0("docker push ", docker_repo))
})
