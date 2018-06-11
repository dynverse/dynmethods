#' Inferring a trajectory inference using [gpfates](https://doi.org/10.1126/sciimmunol.aal2192)
#' 
#' Will generate a trajectory using [gpfates](https://doi.org/10.1126/sciimmunol.aal2192). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/gpfates).
#' 
#' The original code of this method is available [here](https://github.com/Teichlab/GPfates).
#' 
#' The method is described in: [Lönnberg, T., Svensson, V., James, K.R., Fernandez-Ruiz, D., Sebina, I., Montandon, R., Soon, M.S.F., Fogg, L.G., Nair, A.S., Liligeto, U.N., Stubbington, M.J.T., Ly, L.-H., Bagger, F.O., Zwiessele, M., Lawrence, N.D., Souza-Fonseca-Guimaraes, F., Bunn, P.T., Engwerda, C.R., Heath, W.R., Billker, O., Stegle, O., Haque, A., Teichmann, S.A., 2017. Single-cell RNA-seq and computational analysis using temporal mixture modeling resolves TH1/TFHfate bifurcation in malaria. Science Immunology 2, eaal2192.](https://doi.org/10.1126/sciimmunol.aal2192)
#' 
#' @param log_expression_cutoff The log expression cutoff \cr 
#'     numeric; default: 0.5; possible values between 0.5 and 5
#' @param min_cells_expression_cutoff The min expression cutoff \cr 
#'     numeric; default: 0L; possible values between 0 and 20
#' @param ndim Number of dimensions for dimensionality reduction \cr 
#'     integer; default: 2L; possible values between 1 and 5
#' 
#' @return The trajectory model
#' @export
ti_gpfates <- function(
    log_expression_cutoff = 0.5,
    min_cells_expression_cutoff = 0L,
    ndim = 2L
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/gpfates')
  do.call(method, args)
}