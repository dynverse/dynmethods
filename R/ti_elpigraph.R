#' Inferring a trajectory inference using [elpigraph](https://doi.org/https://github.com/Albluca/ElPiGraph.R)
#' 
#' Will generate a trajectory using [elpigraph](https://doi.org/https://github.com/Albluca/ElPiGraph.R). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/elpigraph).
#' 
#' The original code of this method is available [here](https://github.com/Albluca/ElPiGraph.R).
#' 
#' 
#' 
#' @param topology The kind of topology to detect \cr 
#'     discrete; default: "tree"; possible values: cycle, tree, linear
#' @param NumNodes The number of nodes of the principal graph \cr 
#'     integer; default: 50L; possible values between 2 and 1000
#' @param NumEdges The maximum number of edges \cr 
#'     integer; default: 100000L; possible values between 2 and 100000
#' @param InitNodes Number of points to include in the initial graph \cr 
#'     integer; default: 2L; possible values between 2 and 1000
#' @param Mu Controls the elastic energy \cr 
#'     numeric; default: 0.1; possible values between 0.001 and 1
#' @param Lambda Controls the elastic energy \cr 
#'     numeric; default: 0.01; possible values between 0.001 and 1
#' @param MaxNumberOfIterations Maximum number of steps to embed the nodes \cr 
#'     integer; default: 10L; possible values between 1 and 1000
#' @param eps Minimal relative change in the position of the nodes to stop embedment \cr 
#'     numeric; default: 0.01; possible values between 0.001 and 1
#' @param CenterData Should data and initial node positions be centered? \cr 
#' 
#' @return The trajectory model
#' @export
ti_elpigraph <- function(
    topology = "tree",
    NumNodes = 50L,
    NumEdges = 100000L,
    InitNodes = 2L,
    Mu = 0.1,
    Lambda = 0.01,
    MaxNumberOfIterations = 10L,
    eps = 0.01,
    CenterData = FALSE
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/elpigraph')
  do.call(method, args)
}