#' Inferring a trajectory inference using [scimitar](https://doi.org/10.1142/9789813207813_0053)
#' 
#' Will generate a trajectory using [scimitar](https://doi.org/10.1142/9789813207813_0053). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/scimitar).
#' 
#' The original code of this method is available [here](https://github.com/dimenwarper/scimitar).
#' 
#' The method is described in: [CORDERO, P., STUART, J.M., 2016. TRACING CO-REGULATORY NETWORK DYNAMICS IN NOISY, SINGLE-CELL TRANSCRIPTOME TRAJECTORIES. Biocomputing 2017.](https://doi.org/10.1142/9789813207813_0053)
#' 
#' @param covariance_type  \cr 
#'     discrete; default: "diag"; possible values: diag, spherical, full
#' @param degree  \cr 
#'     integer; default: 3L; possible values between 1 and 20
#' @param step_size  \cr 
#'     numeric; default: 0.07; possible values between 0.01 and 0.1
#' @param cov_estimator  \cr 
#'     discrete; default: "corpcor"; possible values: identity, diag, sample, global, glasso, corpcor, average
#' @param cov_reg  \cr 
#'     numeric; default: 0.05; possible values between 0.01 and 0.1
#' @param max_iter  \cr 
#'     integer; default: 3L; possible values between 1 and 20
#' 
#' @return The trajectory model
#' @export
ti_scimitar <- function(
    covariance_type = "diag",
    degree = 3L,
    step_size = 0.07,
    cov_estimator = "corpcor",
    cov_reg = 0.05,
    max_iter = 3L
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/scimitar')
  do.call(method, args)
}