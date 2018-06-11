#' Inferring a trajectory inference using [topslam](https://doi.org/10.1101/057778)
#' 
#' Will generate a trajectory using [topslam](https://doi.org/10.1101/057778). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/topslam).
#' 
#' The original code of this method is available [here](https://github.com/mzwiessele/topslam).
#' 
#' The method is described in: [Zwiessele, M., Lawrence, N.D., 2016. Topslam: Waddington Landscape Recovery for Single Cell Experiments.](https://doi.org/10.1101/057778)
#' 
#' @param n_components The number of components \cr 
#'     integer; default: 2L; possible values between 2 and 10
#' @param n_neighbors The number of neighbors \cr 
#'     integer; default: 10L; possible values between 2 and 100
#' @param linear_dims  \cr 
#'     integer; default: 0L; possible values between 0 and 5
#' @param max_iters The number of iterations to optimize over \cr 
#'     integer; default: 1000L; possible values between 10 and 10000
#' @param dimreds Which dimensionality reductions to use; tSNE, PCA, Spectral, Isomap and/or ICA \cr 
#' 
#' @return The trajectory model
#' @export
ti_topslam <- function(
    n_components = 2L,
    n_neighbors = 10L,
    linear_dims = 0L,
    max_iters = 1000L,
    dimreds = c(TRUE, TRUE, TRUE, TRUE, TRUE)
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/topslam')
  do.call(method, args)
}