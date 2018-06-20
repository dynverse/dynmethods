create_ti_method_chooser <- function(method_function, docker_container) {
  # create arguments
  args <- formals(method_function)
  arg_ids <- names(args)

  # create function
  func <- function(docker = dynwrap::test_docker_installation()) {
    if(docker) {
      purrr::invoke(create_docker_ti_method(docker_container), as.list(environment())[arg_ids])
    } else {
      purrr::invoke(method_function, as.list(environment())[arg_ids])
    }
  }
  formals(func) <- c(formals(func), args)

  func
}
