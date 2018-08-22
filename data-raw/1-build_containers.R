library(tidyverse)
library(furrr)
library(dynmethods)
library(dynwrap)

plan(multiprocess)

definition_files <- list.files("containers", pattern = "definition.yml", recursive = TRUE, full.names = TRUE)
log_dir <- "logs/"
unlink(log_dir, recursive = TRUE)
dir.create(log_dir)

# rebuild all dockers in the 'containers' folder
out <- future_map(definition_files, function(file) {
  # read definition
  definition <- yaml::read_yaml(file)
  docker_repo <- definition$docker_repository
  folder <- stringr::str_replace(file, "/definition.yml$", "")

  # processx::run("docker", args = c("pull", docker_repo), echo = TRUE)

  message("Building ", definition$id)
  err_fun <- function(x, proc) readr::write_lines(x, paste0(log_dir, definition$id, "_stderr.txt"), append = TRUE)
  out_fun <- function(x, proc) readr::write_lines(x, paste0(log_dir, definition$id, "_stdout.txt"), append = TRUE)
  processx::run("docker", args = c("build", folder, "-t", docker_repo), echo = FALSE, stdout_callback = out_fun, stderr_callback = err_fun)
  message("Finished ", definition$id)
})


# rebuild all dockers in the 'containers' folder
out <- map(definition_files, function(file) {
  # read definition
  definition <- yaml::read_yaml(file)
  docker_repo <- definition$docker_repository
  folder <- stringr::str_replace(file, "/definition.yml$", "")

  tryCatch({
    # read example data and possible parameters
    message("Generating example ", definition$id)
    source(paste0(folder, "/example.R"))

    # run method on example, with possible parameters
    message("Running ", definition$id)
    traj <- dynwrap::infer_trajectory(data, docker_repo, parameters = params, verbose = FALSE)

    # if traj is indeed a trajectory, push the docker to dockerhub
    if (dynwrap::is_wrapper_with_trajectory(traj)) {
      message("Pushing ", definition$id)
      processx::run("docker", args = c("push", docker_repo), echo = FALSE)

      data_frame(id = definition$id, built = TRUE, example = TRUE)
    } else {
      message("METHOD ERROR ", definition$id)
      stop("No trajectory found")
    }
  }, error = function(e) {
    message("RUN ERROR ", definition$id, ": ", e$msg)
    data_frame(id = definition$id, built = TRUE, example = FALSE)
  })
})

