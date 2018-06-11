#' Inferring a trajectory inference using [matcher](https://doi.org/10.1186/s13059-017-1269-0)
#' 
#' Will generate a trajectory using [matcher](https://doi.org/10.1186/s13059-017-1269-0). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/matcher).
#' 
#' The original code of this method is available [here](https://github.com/jw156605/MATCHER).
#' 
#' The method is described in: [Welch, J.D., Hartemink, A.J., Prins, J.F., 2017. MATCHER: manifold alignment reveals correspondence between single cell transcriptome and epigenome dynamics. Genome Biology 18.](https://doi.org/10.1186/s13059-017-1269-0)
#' 
#' @param quantiles Quantiles How many quantiles to use when computing warp functions (integer) \cr 
#'     integer; default: 50L; possible values between 2 and 500
#' @param method Gaussian process regression or linear interpolation? ("gp" or "linear) \cr 
#'     discrete; default: "linear"; possible values: gp, linear
#' 
#' @return The trajectory model
#' @export
ti_matcher <- function(
    quantiles = 50L,
    method = "linear"
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/matcher')
  do.call(method, args)
}