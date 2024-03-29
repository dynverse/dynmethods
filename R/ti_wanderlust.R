######################################### DO NOT EDIT! #########################################
#### This file is automatically generated from data-raw/2-generate_r_code_from_containers.R ####
################################################################################################

#' @title Wanderlust
#' 
#' @description
#' Will generate a trajectory using
#' [Wanderlust](https://doi.org/10.1016/j.cell.2014.04.005).
#' 
#' This method was wrapped inside a
#' [container](https://github.com/dynverse/ti_wanderlust).
#' The original code of this method is available
#' [here](https://github.com/ManuSetty/wishbone).
#' 
#' @references Bendall, S.C., Davis, K.L., Amir, E.D., Tadmor, M.D., Simonds,
#' E.F., Chen, T.J., Shenfeld, D.K., Nolan, G.P., Pe’er, D., 2014. Single-Cell
#' Trajectory Detection Uncovers Progression and Regulatory Coordination in Human
#' B Cell Development. Cell 157, 714–725.
#' 
#' @param normalise . Default: TRUE. Format: logical.
#' @param knn K-nearest neighbours for diffusion. Domain: U(15, 100). Default: 25.
#' Format: integer.
#' @param n_diffusion_components Number of diffusion components. Domain: U(3, 20).
#' Default: 3. Format: integer.
#' @param n_pca_components Number of pca components. Domain: U(15, 100). Default:
#' 30. Format: integer.
#' @param k K parameter. Domain: U(15, 100). Default: 25. Format: integer.
#' @param num_waypoints Number of waypoints. Domain: U(100, 500). Default: 250.
#' Format: integer.
#' @param epsilon Epsilon. Domain: U(0.1, 5). Default: 1. Format: numeric.
#' 
#' @keywords method
#' 
#' @return A TI method wrapper to be used together with
#' \code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
#' @export
ti_wanderlust <- function(
    normalise = TRUE,
    knn = 25L,
    n_diffusion_components = 3L,
    n_pca_components = 30L,
    k = 25L,
    num_waypoints = 250L,
    epsilon = 1L
) {
  method_choose_backend(
    package_repository = NULL,
    package_name = NULL,
    function_name = NULL,
    package_version = NULL,
    container_id = "dynverse/ti_wanderlust:v0.9.9.01"
  )(
    normalise = normalise,
    knn = knn,
    n_diffusion_components = n_diffusion_components,
    n_pca_components = n_pca_components,
    k = k,
    num_waypoints = num_waypoints,
    epsilon = epsilon
  )
}

