#' Inferring a trajectory inference using [grandprix](https://doi.org/10.1101/227843)
#' 
#' Will generate a trajectory using [grandprix](https://doi.org/10.1101/227843). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/grandprix).
#' 
#' The original code of this method is available [here](https://github.com/ManchesterBioinference/GrandPrix).
#' 
#' The method is described in: [Ahmed, S., Rattray, M., Boukouvalas, A., 2017. GrandPrix: Scaling up the Bayesian GPLVM for single-cell data.](https://doi.org/10.1101/227843)
#' 
#' @param n_latent_dims  \cr 
#'     integer; default: 2L; possible values between 1 and 10
#' @param n_inducing_points  \cr 
#'     integer; default: 40L; possible values between 10 and 500
#' @param latent_prior_var  \cr 
#' @param latent_var  \cr 
#' 
#' @return The trajectory model
#' @export
ti_grandprix <- function(
    n_latent_dims = 2L,
    n_inducing_points = 40L,
    latent_prior_var = 0.1,
    latent_var = 0.028
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/grandprix')
  do.call(method, args)
}