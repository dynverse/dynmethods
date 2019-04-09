
#' @importFrom dynwrap create_ti_method_container
method_choose_backend <- function(r_remote, r_function, container_id) {
  correct_backends <- c("r_wrapper", "container")

  # first choose backend: either with an option, using a selection menu (if interactive) or defaults to r_wrapper
  if (!is.null(getOption("dynwrap_backend"))) {
    backend <- getOption("dynwrap_backend")
    if (!backend %in% correct_backends) {
      stop("Invalid dynwrap_backend option: ", backend)
    }
  } else if (interactive()) {

    # prompt user for container vs r_wrapper
    title <- paste(
      paste0("You can run this method as an ", crayon::bold("R wrapper"), " (1, default) or as a ", crayon::bold("container"), " (2)"),
      "Which do you want to use? This option will be saved in options(dynwrap_backend = c('r_wrapper', 'container'))",
      "1: R wrapper [default]",
      "2: Container",
      sep = "\n"
    )

    message(title)
    answer <- readline("")

    if (answer == "2") {
      backend <- "container"
    } else {
      backend <- "r_wrapper"
    }

    # set option
    options(dynwrap_backend = backend)
  } else {
    message("dynwrap_backend option not set, assuming that you want to run this method as an R wrapper")

    backend <- "r_wrapper"
  }

  # now return method
  if (backend == "r_wrapper") {
    install_github_tagged_version(r_remote)
    package <- parse_github_repo_spec(r_remote)$package
    getFromNamespace(r_function, package)

  } else if (backend == "container") {
    dynwrap::create_ti_method_container(container_id = container_id)
  }
}
