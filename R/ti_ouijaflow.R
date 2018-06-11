#' Inferring a trajectory inference using [ouijaflow](https://doi.org/10.1101/060442)
#' 
#' Will generate a trajectory using [ouijaflow](https://doi.org/10.1101/060442). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/ouijaflow).
#' 
#' The original code of this method is available [here](https://github.com/kieranrcampbell/ouija).
#' 
#' The method is described in: [Campbell, K.R., Yau, C., 2016. A descriptive marker gene approach to single-cell pseudotime inference.](https://doi.org/10.1101/060442)
#' 
#' @param iter  \cr 
#'     integer; default: 1000L; possible values between 2 and 50000
#' 
#' @return The trajectory model
#' @export
ti_ouijaflow <- function(
    iter = 1000L
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/ouijaflow')
  do.call(method, args)
}