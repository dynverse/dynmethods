# this should probably be moved to dynwrap
# and be merged with other functions like
# * create_image_ti_method
# * create_docker_ti_method
# * create_singularity_ti_method
# * create_ti_method

#' Creating a TI method from a docker repository
#'
#' Supports both docker and singularity as a backend.
#'
#' @param docker_repository The id of the docker repository.
#' @param run_environment In which environment to run the method, can be `"docker"` or `"singularity"`.
#' @param ... Default parameters for the method.
create_container_ti_method <- function(docker_repository, run_environment = NULL, ...) {
  if (is.null(run_environment)) {
    run_environment <- getOption("dynwrap_run_environment")
  }
  if (is.null(run_environment)) {
    run_environment <- "docker"
  }

  # choose environments
  if (run_environment == "docker") {
    method <- create_docker_ti_method(docker_repository)
  } else if (run_environment == 'singularity') {
    method <- create_singularity_ti_method(docker_repository, ".simg")
  }

  method(...)
}
