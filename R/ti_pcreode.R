#' Inferring a trajectory inference using [pCreode](https://doi.org/10.1016/j.cels.2017.10.012)
#' 
#' Will generate a trajectory using [pCreode](https://doi.org/10.1016/j.cels.2017.10.012). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/pcreode).
#' 
#' The original code of this method is available [here](https://github.com/KenLauLab/pCreode).
#' 
#' The method is described in: [Herring, C.A., Banerjee, A., McKinley, E.T., Simmons, A.J., Ping, J., Roland, J.T., Franklin, J.L., Liu, Q., Gerdes, M.J., Coffey, R.J., Lau, K.S., 2018. Unsupervised Trajectory Analysis of Single-Cell RNA-Seq and Imaging Data Reveals Alternative Tuft Cell Origins in the Gut. Cell Systems 6, 37–51.e9.](https://doi.org/10.1016/j.cels.2017.10.012)
#' 
#' @param n_pca_components  \cr 
#'     integer; default: 3L; possible values between 2 and 10
#' @param radius  \cr 
#'     numeric; default: 1L; possible values between 0.01 and 10
#' @param noise  \cr 
#'     numeric; default: 8L; possible values between 1 and 20
#' @param target  \cr 
#'     numeric; default: 25L; possible values between 5 and 100
#' @param num_runs  \cr 
#'     integer; default: 10L; possible values between 10 and 1000
#' 
#' @return The trajectory model
#' @export
ti_pcreode <- function(
    n_pca_components = 3L,
    radius = 1L,
    noise = 8L,
    target = 25L,
    num_runs = 10L
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/pcreode')
  do.call(method, args)
}