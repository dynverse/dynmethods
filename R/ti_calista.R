######################################### DO NOT EDIT! #########################################
#### This file is automatically generated from data-raw/2-generate_r_code_from_containers.R ####
################################################################################################

#' @title CALISTA
#' 
#' @description
#' Will generate a trajectory using [CALISTA](https://doi.org/10.1101/257550).
#' 
#' This method was wrapped inside a
#' [container](https://github.com/dynverse/ti_calista).
#' The original code of this method is available
#' [here](https://github.com/CABSEL/CALISTA).
#' 
#' @references Gao, N.P., Hartmann, T., Fang, T., Gunawan, R., 2018. CALISTA:
#' Clustering and Lineage Inference in Single-Cell Transcriptional Analysis.
#' 
#' @param runs Number of independent runs of greedy algorithm. Domain: U(20, 100).
#' Default: 50. Format: integer.
#' @param max_iter Number of iterations in greedy algorithm. Domain: U(20, 400).
#' Default: 100. Format: integer.
#' 
#' @keywords method
#' 
#' @return A TI method wrapper to be used together with
#' \code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
#' @export
ti_calista <- function(
    runs = 50L,
    max_iter = 100L
) {
  create_ti_method_container(container_id = "dynverse/ti_calista:v0.9.9.01")(
    runs = runs,
    max_iter = max_iter
  )
}

