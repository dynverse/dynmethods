######################################### DO NOT EDIT! #########################################
#### This file is automatically generated from data-raw/2-generate_r_code_from_containers.R ####
################################################################################################

#' @title Identity
#' 
#' @description
#' Will generate a trajectory using Identity.
#' 
#' This method was wrapped inside a
#' [container](https://github.com/dynverse/ti_identity).
#' 
#' 
#' 
#' 
#' 
#' @keywords method
#' 
#' @return A TI method wrapper to be used together with
#' \code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
#' @export
ti_identity <- function(
    
) {
  method_choose_backend(
    package_repository = NULL,
    package_name = NULL,
    function_name = NULL,
    package_version = NULL,
    container_id = "dynverse/ti_identity:v0.9.9.01"
  )(
    
  )
}

