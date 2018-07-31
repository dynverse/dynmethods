library(tidyverse)
library(furrr)

plan(multiprocess)

definition_files <- list.files("containers", pattern = "definition.yml", recursive = TRUE, full.names = TRUE)

# rebuild all dockers in the 'containers' folder
out <- future_map(definition_files, function(file) {
  definition <- yaml::read_yaml(file)
  docker_repo <- definition$docker_repository
  folder <- stringr::str_replace(file, "/definition.yml$", "")

  # processx::run("docker", args = c("pull", docker_repo), echo = TRUE)
  processx::run("docker", args = c("build", folder, "-t", docker_repo), echo = TRUE)

  tryCatch({
    rm(data, params)
    source(stringr::str_replace(file, "/definition.yml$", "/example.R"))

    traj <- dynwrap::infer_trajectory(data, definition$id, parameters = {if (exists("params")) params else NULL}, verbose = FALSE)

    if (dynwrap::is_wrapper_with_trajectory(traj)) {
      processx::run("docker", args = c("push", docker_repo), echo = TRUE)

      data_frame(id = definition$id, built = TRUE, example = TRUE)
    } else {
      stop("No trajectory found")
    }
  }, error = function(e) {
    data_frame(id = definition$id, built = TRUE, example = FALSE)
  })
})
