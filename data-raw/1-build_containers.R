library(tidyverse)
library(furrr)
library(dynmethods)
library(dynwrap)

plan(multiprocess)

files <- list.files("containers", pattern = "definition.yml", recursive = TRUE, full.names = TRUE)
log_dir <- "logs/"
unlink(log_dir, recursive = TRUE)
dir.create(log_dir)

# rebuild all dockers in the 'containers' folder
out <- map(files, function(file) {
  folder <- stringr::str_replace(file, "/definition.yml$", "")

  definition <- yaml::read_yaml(file)
  method_id <- definition$id
  docker_repo <- definition$docker_repository

  err_fun <- function(x, proc) readr::write_lines(x, paste0(log_dir, method_id, "_stderr.txt"), append = TRUE)
  out_fun <- function(x, proc) readr::write_lines(x, paste0(log_dir, method_id, "_stdout.txt"), append = TRUE)

  message(method_id, ": Pull remote docker")
  processx::run("docker", args = c("pull", docker_repo), echo = FALSE, stdout_callback = out_fun, stderr_callback = err_fun)
  digest <- tryCatch(dynwrap:::.container_get_digests(docker_repo, "docker")$digest, error = function(e) NA)

  message(method_id, ": Building new docker")
  processx::run("docker", args = c("build", folder, "-t", docker_repo), echo = FALSE, stdout_callback = out_fun, stderr_callback = err_fun)

  new_digest <- tryCatch(dynwrap:::.container_get_digests(docker_repo, "docker")$digest, error = function(e) NA)

  if (!identical(new_digest, digest)) {
    tryCatch({
      message(method_id, ": New version found; testing method")

      # read example data and possible parameters
      message(method_id, ": Generating example and params")
      source(paste0(folder, "/example.R"))

      message(method_id, ": Testing par_set")
      method <- dynwrap::create_ti_method_with_container(docker_repo)()
      par_set <- method$par_set
      design <- ParamHelpers::generateDesign(3, par_set, trafo = TRUE)
      design <- ParamHelpers::generateDesignOfDefaults(par_set, trafo = TRUE)

      # run method on example, with possible parameters
      message(method_id, ": Executing method")
      if (method_id != "error") {
        traj <- dynwrap::infer_trajectory(data, docker_repo, parameters = params, verbose = FALSE)
      }

      # if traj is indeed a trajectory, push the docker to dockerhub
      message(method_id, ": OK! Pushing to docker")
      processx::run("docker", args = c("push", docker_repo), echo = FALSE)

    }, error = function(e) {
      data_frame(id = method_id, built = TRUE, example = FALSE)
    })
  } else {
    message(method_id, ": No changes detected")
  }

  message(method_id, ": Finished ")
})
