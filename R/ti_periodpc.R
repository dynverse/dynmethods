######################################### DO NOT EDIT! #########################################
#### This file is automatically generated from data-raw/2-generate_r_code_from_containers.R ####
################################################################################################

#' @title Periodic PrinCurve
#' 
#' @description
#' Will generate a trajectory using Periodic PrinCurve.
#' 
#' This method was wrapped inside a
#' [container](https://github.com/dynverse/ti_periodpc).
#' 
#' 
#' 
#' 
#' @param ndim . Domain: U(2, 10). Default: 3. Format: integer.
#' @param maxit . Domain: U(0, 100). Default: 10. Format: integer.
#' 
#' @keywords method
#' 
#' @return A TI method wrapper to be used together with
#' \code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
#' @export
ti_periodpc <- function(
    ndim = 3L,
    maxit = 10L
) {
  method_choose_backend(
    package_repository = NULL,
    package_name = NULL,
    function_name = NULL,
    package_version = NULL,
    container_id = "dynverse/ti_periodpc:v0.9.9.01"
  )(
    ndim = ndim,
    maxit = maxit
  )
}

