#' Methods for working with the dynmethods docker
#'
#' Make sure docker is installed and available through path, see \url{https://docs.docker.com/install}.
#' After running `start_dynmethods_docker`, methods can be run on the docker using \code{x \%<-\% ...}
#'
#' @examples
#' start_dynmethods_docker()
#' tasks <- dyntoy::toy_tasks[1, ]
#' models %<-% infer_trajectory(tasks, description_compone())
#'
#' @export
start_dynmethods_docker <- function() {
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

#' @importFrom stringr str_glue
#' @rdname start_dynmethods_docker
stop_dynmethods_docker <- function() {
  to_stop <- system("docker ps -a | grep dynmethods | awk '{print $1}'", intern=T)

  if (length(to_stop)) {
    system(str_glue("docker stop {to_stop}"))
    system(str_glue("docker rm {to_stop}"))
  }
}

check_dynmethods_docker_running <- function() {
  system("docker ps -a | grep dynmethods | awk '{print $1}'", intern=T)
}

#' @rdname start_dynmethods_docker
start_dynmethods_docker <- function() {
  if(length(check_dynmethods_docker_running())) {
    stop_dynmethods_docker()
  }
  cl <- future::makeClusterPSOCK(
    "localhost",
    rscript = c(
      "docker", "run", "--network=host", "--name=dynmethods", "dynverse/dynmethods",
      "Rscript"
    ),
    rscript_args = c(
      "-e",
      shQuote("library(dynmethods)")
    ),
    connectTimeout=10
  )
  future::plan(future::cluster, workers = cl)
}

#' @importFrom future %<-%
#' @export
future::`%<-%`

