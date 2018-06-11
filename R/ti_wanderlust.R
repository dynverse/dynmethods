#' Inferring a trajectory inference using [wanderlust](https://doi.org/10.1016/j.cell.2014.04.005)
#' 
#' Will generate a trajectory using [wanderlust](https://doi.org/10.1016/j.cell.2014.04.005). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/wanderlust).
#' 
#' The original code of this method is available [here](https://www.c2b2.columbia.edu/danapeerlab/html/wanderlust.html).
#' 
#' The method is described in: [Bendall, S.C., Davis, K.L., Amir, E.D., Tadmor, M.D., Simonds, E.F., Chen, T.J., Shenfeld, D.K., Nolan, G.P., Pe’er, D., 2014. Single-Cell Trajectory Detection Uncovers Progression and Regulatory Coordination in Human B Cell Development. Cell 157, 714–725.](https://doi.org/10.1016/j.cell.2014.04.005)
#' 
#' @param branch Whether to allow a single bifurcation within the trajectory (wishbone versus wanderlust) \cr 
#' @param epsilon Epsilon \cr 
#'     numeric; default: 1L; possible values between 0.1 and 10
#' @param k K parameter \cr 
#'     integer; default: 15L; possible values between 2 and 100
#' @param knn K-nearest neighbours for diffusion \cr 
#'     integer; default: 15L; possible values between 2 and 100
#' @param n_diffusion_components Number of diffusion components \cr 
#'     integer; default: 2L; possible values between 2 and 20
#' @param n_pca_components Number of pca components \cr 
#'     integer; default: 15L; possible values between 2 and 30
#' @param normalise  \cr 
#' @param num_waypoints Number of waypoints \cr 
#'     integer; default: 250L; possible values between 2 and 500
#' 
#' @return The trajectory model
#' @export
ti_wanderlust <- function(
    branch = FALSE,
    epsilon = 1L,
    k = 15L,
    knn = 15L,
    n_diffusion_components = 2L,
    n_pca_components = 15L,
    normalise = TRUE,
    num_waypoints = 250L
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/wanderlust')
  do.call(method, args)
}