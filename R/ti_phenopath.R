######################################### DO NOT EDIT! #########################################
#### This file is automatically generated from data-raw/2-generate_r_code_from_containers.R ####
################################################################################################

#' @title PhenoPath
#' 
#' @description
#' Will generate a trajectory using [PhenoPath](https://doi.org/10.1101/159913).
#' 
#' This method was wrapped inside a
#' [container](https://github.com/dynverse/ti_phenopath).
#' The original code of this method is available
#' [here](https://github.com/kieranrcampbell/phenopath).
#' 
#' @references Campbell, K., Yau, C., 2017. Uncovering genomic trajectories with
#' heterogeneous genetic and environmental backgrounds across single-cells and
#' populations.
#' 
#' @param thin The number of iterations to wait each time before re-calculating
#' the elbo. Domain: U(2, 500). Default: 40. Format: integer.
#' @param z_init The initialisation of the latent trajectory. Should be one of* A
#' positive integer describing which principal component of the data should be
#' used for initialisation (default 1), or * The text character `"random"`, for
#' random initialisation from a standard normal distribution. Domain: {1, 2, 3, 4,
#' 5, random}. Default: 1. Format: character.
#' @param model_mu Logical - should a gene-specific intercept term be modelled?.
#' Default: FALSE. Format: logical.
#' @param scale_y Logical - should the expression matrix be centre scaled?.
#' Default: TRUE. Format: logical.
#' 
#' @keywords method
#' 
#' @return A TI method wrapper to be used together with
#' \code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
#' @export
ti_phenopath <- function(
    thin = 40L,
    z_init = "1",
    model_mu = FALSE,
    scale_y = TRUE
) {
  method_choose_backend(
    package_repository = NULL,
    package_name = NULL,
    function_name = NULL,
    package_version = NULL,
    container_id = "dynverse/ti_phenopath:v0.9.9.01"
  )(
    thin = thin,
    z_init = z_init,
    model_mu = model_mu,
    scale_y = scale_y
  )
}

