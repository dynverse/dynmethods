#' Inferring a trajectory inference using [wishbone](https://doi.org/10.1038/nbt.3569)
#' 
#' Will generate a trajectory using [wishbone](https://doi.org/10.1038/nbt.3569). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/wishbone).
#' 
#' The original code of this method is available [here](https://github.com/ManuSetty/wishbone).
#' 
#' The method is described in: [Setty, M., Tadmor, M.D., Reich-Zeliger, S., Angel, O., Salame, T.M., Kathail, P., Choi, K., Bendall, S., Friedman, N., Pe’er, D., 2016. Wishbone identifies bifurcating developmental trajectories from single-cell data. Nature Biotechnology 34, 637–645.](https://doi.org/10.1038/nbt.3569)
#' 
#' @param normalise  \cr 
#' @param knn K-nearest neighbours for diffusion \cr 
#'     integer; default: 15L; possible values between 2 and 100
#' @param n_diffusion_components Number of diffusion components \cr 
#'     integer; default: 2L; possible values between 2 and 20
#' @param n_pca_components Number of pca components \cr 
#'     integer; default: 15L; possible values between 2 and 30
#' @param k K parameter \cr 
#'     integer; default: 15L; possible values between 2 and 100
#' @param num_waypoints Number of waypoints \cr 
#'     integer; default: 250L; possible values between 2 and 500
#' @param epsilon Epsilon \cr 
#'     numeric; default: 1L; possible values between 0.1 and 10
#' @param branch Whether to allow a single bifurcation within the trajectory (wishbone versus wanderlust) \cr 
#' 
#' @return The trajectory model
#' @export
ti_wishbone <- function(
    normalise = TRUE,
    knn = 15L,
    n_diffusion_components = 2L,
    n_pca_components = 15L,
    k = 15L,
    num_waypoints = 250L,
    epsilon = 1L,
    branch = TRUE
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/wishbone')
  do.call(method, args)
}