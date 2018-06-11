#' Inferring a trajectory inference using [scuba](https://doi.org/10.1073/pnas.1408993111)
#' 
#' Will generate a trajectory using [scuba](https://doi.org/10.1073/pnas.1408993111). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/scuba).
#' 
#' The original code of this method is available [here](https://github.com/gcyuan/SCUBA).
#' 
#' The method is described in: [Marco, E., Karp, R.L., Guo, G., Robson, P., Hart, A.H., Trippa, L., Yuan, G.-C., 2014. Bifurcation analysis of single-cell gene expression data reveals epigenetic landscape. Proceedings of the National Academy of Sciences 111, E5643–E5650.](https://doi.org/10.1073/pnas.1408993111)
#' 
#' @param rigorous_gap_stats Whether to use rigorous gap statistics to determine number of clusters \cr 
#' @param N_dim Number of TSNE dimensions \cr 
#'     integer; default: 2L; possible values between 2 and 3
#' @param low_gene_threshold Threshold value for genes of low expression levels \cr 
#'     numeric; default: 1L; possible values between 0 and 5
#' @param low_gene_fraction_max Maximum fraction of lowly-expressed cells allowed for each gene \cr 
#'     numeric; default: 0.7; possible values between 0 and 1
#' @param min_split Lower threshold on the number of cells in a cluster for this cluster to be split. \cr 
#'     integer; default: 15L; possible values between 1 and 100
#' @param min_percentage_split Minimum fraction of cells in the smaller cluster during a bifurcation. \cr 
#'     numeric; default: 0.25; possible values between 0 and 1
#' 
#' @return The trajectory model
#' @export
ti_scuba <- function(
    rigorous_gap_stats = TRUE,
    N_dim = 2L,
    low_gene_threshold = 1L,
    low_gene_fraction_max = 0.7,
    min_split = 15L,
    min_percentage_split = 0.25
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/scuba')
  do.call(method, args)
}