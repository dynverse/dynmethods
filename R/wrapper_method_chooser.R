create_ti_method_chooser <- function(method_function, docker_container) {
  # create arguments
  args <- formals(method_function)
  arg_ids <- names(args)

  # create function
  func <- function(
    run_environment = ifelse(is.null(getOption("dynwrap_run_environment")), "local", getOption("dynwrap_run_environment"))
  ) {
    # choose environments
    if (run_environment == "docker") {
      purrr::invoke(create_docker_ti_method(docker_container), as.list(environment())[arg_ids])
    } else if (run_environment == "singularity") {
      purrr::invoke(create_singularity_ti_method(paste0(docker_container, ".simg")), as.list(environment())[arg_ids])
    } else {
      purrr::invoke(method_function, as.list(environment())[arg_ids])
    }
  }
  formals(func) <- c(formals(func), args)

  func
}
