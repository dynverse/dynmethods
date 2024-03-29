######################################### DO NOT EDIT! #########################################
#### This file is automatically generated from data-raw/2-generate_r_code_from_containers.R ####
################################################################################################

#' @title MFA
#' 
#' @description
#' Will generate a trajectory using
#' [MFA](https://doi.org/10.12688/wellcomeopenres.11087.1).
#' 
#' This method was wrapped inside a
#' [container](https://github.com/dynverse/ti_mfa).
#' The original code of this method is available
#' [here](https://github.com/kieranrcampbell/mfa).
#' 
#' @references Campbell, K.R., Yau, C., 2017. Probabilistic modeling of
#' bifurcations in single-cell gene expression data using a Bayesian mixture of
#' factor analyzers. Wellcome Open Research 2, 19.
#' 
#' @param iter Number of MCMC iterations. Domain: U(20, 5000). Default: 2000.
#' Format: integer.
#' @param thin MCMC samples to thin. Domain: U(1, 20). Default: 1. Format:
#' integer.
#' @param pc_initialise Which principal component to initialise pseudotimes to.
#' Domain: U(1, 5). Default: 1. Format: integer.
#' @param prop_collapse Proportion of Gibbs samples which should marginalise over
#' c. Domain: U(0, 1). Default: 0. Format: numeric.
#' @param scale_input Logical. If true, input is scaled to have mean 0 variance 1.
#' Default: TRUE. Format: logical.
#' @param zero_inflation Logical, should zero inflation be enabled?. Default:
#' FALSE. Format: logical.
#' 
#' @keywords method
#' 
#' @return A TI method wrapper to be used together with
#' \code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
#' @export
ti_mfa <- function(
    iter = 2000L,
    thin = 1L,
    pc_initialise = 1L,
    prop_collapse = 0L,
    scale_input = TRUE,
    zero_inflation = FALSE
) {
  method_choose_backend(
    package_repository = NULL,
    package_name = NULL,
    function_name = NULL,
    package_version = NULL,
    container_id = "dynverse/ti_mfa:v0.9.9.01"
  )(
    iter = iter,
    thin = thin,
    pc_initialise = pc_initialise,
    prop_collapse = prop_collapse,
    scale_input = scale_input,
    zero_inflation = zero_inflation
  )
}

