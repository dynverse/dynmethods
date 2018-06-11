#' Inferring a trajectory inference using [paga](https://doi.org/10.1101/208819)
#' 
#' Will generate a trajectory using [paga](https://doi.org/10.1101/208819). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/paga).
#' 
#' The original code of this method is available [here](https://github.com/theislab/graph_abstraction).
#' 
#' The method is described in: [Wolf, F.A., Hamey, F., Plass, M., Solana, J., Dahlin, J.S., Gottgens, B., Rajewsky, N., Simon, L., Theis, F.J., 2017. Graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells.](https://doi.org/10.1101/208819)
#' 
#' @param n_neighbors Number of neighbours for knn \cr 
#'     integer; default: 30L; possible values between 1 and 100
#' @param n_comps Number of principal components \cr 
#'     integer; default: 50L; possible values between 0 and 100
#' @param resolution Resolution of louvain clustering, which determines the granularity of the clustering. Higher values will result in more clusters. \cr 
#'     numeric; default: 2.5; possible values between 0.1 and 10
#' 
#' @return The trajectory model
#' @export
ti_paga <- function(
    n_neighbors = 30L,
    n_comps = 50L,
    resolution = 2.5
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/paga')
  do.call(method, args)
}