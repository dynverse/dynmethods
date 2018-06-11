#' Inferring a trajectory inference using [CellRouter](https://doi.org/10.1038/s41467-018-03214-y)
#' 
#' Will generate a trajectory using [CellRouter](https://doi.org/10.1038/s41467-018-03214-y). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/cellrouter).
#' 
#' The original code of this method is available [here](https://github.com/edroaldo/cellrouter).
#' 
#' The method is described in: [Lummertz da Rocha, E., Rowe, R.G., Lundin, V., Malleshaiah, M., Jha, D.K., Rambo, C.R., Li, H., North, T.E., Collins, J.J., Daley, G.Q., 2018. Reconstruction of complex single-cell trajectories using CellRouter. Nature Communications 9.](https://doi.org/10.1038/s41467-018-03214-y)
#' 
#' @param ndim_pca Number of principlal components to compute \cr 
#'     integer; default: 20L; possible values between 2 and 100
#' @param ndim_tsne Number of tsne dimensions to compute \cr 
#'     integer; default: 11L; possible values between 2 and 100
#' @param max_iter Maximal number of tsne iterations \cr 
#'     integer; default: 1000L; possible values between 100 and 100000
#' @param cluster_method Method to use for clustering \cr 
#'     discrete; default: "graph.clustering"; possible values: graph.clustering, model.clustering
#' @param k_clustering Number of nearest neighbors to build a k-nearest neighbors graph for clustering \cr 
#'     integer; default: 20L; possible values between 2 and 1000
#' @param ndim_pca_clustering Number of PCA dimensions used for k-nearest neighbors graph for clustering \cr 
#'     integer; default: 20L; possible values between 2 and 100
#' @param k_knn Number of nearest neighbors to build a k-nearest neighbors graph for knn \cr 
#'     integer; default: 10L; possible values between 2 and 1000
#' @param ndim_pca_knn Number of PCA dimensions used for knn \cr 
#'     integer; default: 20L; possible values between 2 and 100
#' @param sim_type Similarity type for knn \cr 
#'     discrete; default: "jaccard"; possible values: jaccard
#' @param distance_method_paths Distance method for paths \cr 
#'     discrete; default: "graph"; possible values: euclidean, maximum, manhattan, canberra, binary, graph
#' @param ranks How to rank the paths \cr 
#'     discrete; default: "rank"; possible values: path_cost, path_flow, rank, length
#' @param num_cells Trajectories should contain at least num.cells \cr 
#'     integer; default: 3L; possible values between 3 and 100
#' @param neighs The size of the neighborhood in kNN graph used to smoothen kinetic profiles \cr 
#'     integer; default: 3L; possible values between 2 and 100
#' 
#' @return The trajectory model
#' @export
ti_cellrouter <- function(
    ndim_pca = 20L,
    ndim_tsne = 11L,
    max_iter = 1000L,
    cluster_method = "graph.clustering",
    k_clustering = 20L,
    ndim_pca_clustering = 20L,
    k_knn = 10L,
    ndim_pca_knn = 20L,
    sim_type = "jaccard",
    distance_method_paths = "graph",
    ranks = "rank",
    num_cells = 3L,
    neighs = 3L
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/cellrouter')
  do.call(method, args)
}