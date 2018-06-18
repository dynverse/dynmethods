create_ti_method_chooser <- function(method, docker_container) {
  # create arguments
  args <- formals(method)
  arg_ids <- names(args)

  # create function
  func <- function(docker = "stevedore" %in% rownames(installed.packages())) {
    if(docker) {
      invoke(create_docker_ti_method(docker_container), as.list(environment())[arg_ids])
    } else {
      invoke(method, as.list(environment())[arg_ids])
    }
  }
  formals(func) <- c(formals(func), args)

  func
}
