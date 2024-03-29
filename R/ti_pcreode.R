######################################### DO NOT EDIT! #########################################
#### This file is automatically generated from data-raw/2-generate_r_code_from_containers.R ####
################################################################################################

#' @title pCreode
#' 
#' @description
#' Will generate a trajectory using
#' [pCreode](https://doi.org/10.1016/j.cels.2017.10.012).
#' 
#' This method was wrapped inside a
#' [container](https://github.com/dynverse/ti_pcreode).
#' The original code of this method is available
#' [here](https://github.com/KenLauLab/pCreode).
#' 
#' @references Herring, C.A., Banerjee, A., McKinley, E.T., Simmons, A.J., Ping,
#' J., Roland, J.T., Franklin, J.L., Liu, Q., Gerdes, M.J., Coffey, R.J., Lau,
#' K.S., 2018. Unsupervised Trajectory Analysis of Single-Cell RNA-Seq and Imaging
#' Data Reveals Alternative Tuft Cell Origins in the Gut. Cell Systems 6,
#' 37–51.e9.
#' 
#' @param n_pca_components . Domain: U(2, 10). Default: 3. Format: integer.
#' @param num_runs . Domain: e^U(2.30, 4.61). Default: 10. Format: integer.
#' 
#' @keywords method
#' 
#' @return A TI method wrapper to be used together with
#' \code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
#' @export
ti_pcreode <- function(
    n_pca_components = 3L,
    num_runs = 10L
) {
  method_choose_backend(
    package_repository = NULL,
    package_name = NULL,
    function_name = NULL,
    package_version = NULL,
    container_id = "dynverse/ti_pcreode:v0.9.9.01"
  )(
    n_pca_components = n_pca_components,
    num_runs = num_runs
  )
}

