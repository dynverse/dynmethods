#' Stop the dynmethods docker
#'
#' @export
stop_dynmethods_docker <- function() {
  to_stop <- system("docker ps -a | grep dynmethods | awk '{print $1}'", intern=T)

  if(length(to_stop)) {
    system(str_glue("docker stop {to_stop}"))
    system(str_glue("docker rm {to_stop}"))
  }
}

#' Start the dynmethods docker
#'
#' Make sure docker is installed and available through path
#'
#' @export
start_dynmethods_docker <- function() {
  stop_dynmethods_docker()

  cl <- future::makeClusterPSOCK(
    "localhost",
    rscript = c(
      "docker", "run", "--network=host", "--name=dynmethods", "dynverse/dynmethods",
      "Rscript"
    ),
    rscript_args = c(
      "-e",
      shQuote("library(dynmethods)")
    )
  )
  future::plan(future::cluster, workers = cl)
}

#' @importFrom future %<-%
#' @export
future::`%<-%`
