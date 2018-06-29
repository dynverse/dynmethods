#' Inferring a trajectory inference using Angle
#' 
#' Will generate a trajectory using Angle. This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/angle).
#' 
#' This methods was first wrapped inside R, see [ti_angle]
#' 
#' 
#' 
#' 
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_angle <- create_ti_method_chooser(ti_angle, 'dynverse/angle')




#' Inferring a trajectory inference using [CellRouter](https://doi.org/10.1038/s41467-018-03214-y)
#' 
#' Will generate a trajectory using [CellRouter](https://doi.org/10.1038/s41467-018-03214-y). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/cellrouter).
#' 
#' 
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




#' Inferring a trajectory inference using [cellTree with gibbs](https://doi.org/10.1186/s12859-016-1175-6)
#' 
#' Will generate a trajectory using [cellTree with gibbs](https://doi.org/10.1186/s12859-016-1175-6). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/celltree_gibbs).
#' 
#' This methods was first wrapped inside R, see [ti_celltree_gibbs]
#' 
#' 
#' 
#' The method is described in: [duVerle, D.A., Yotsukura, S., Nomura, S., Aburatani, H., Tsuda, K., 2016. CellTree: an R/bioconductor package to infer the hierarchical structure of cell populations from single-cell RNA-seq data. BMC Bioinformatics 17.](https://doi.org/10.1186/s12859-016-1175-6)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_celltree_gibbs <- create_ti_method_chooser(ti_celltree_gibbs, 'dynverse/celltree_gibbs')




#' Inferring a trajectory inference using [cellTree with maptpx](https://doi.org/10.1186/s12859-016-1175-6)
#' 
#' Will generate a trajectory using [cellTree with maptpx](https://doi.org/10.1186/s12859-016-1175-6). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/celltree_maptpx).
#' 
#' This methods was first wrapped inside R, see [ti_celltree_maptpx]
#' 
#' 
#' 
#' The method is described in: [duVerle, D.A., Yotsukura, S., Nomura, S., Aburatani, H., Tsuda, K., 2016. CellTree: an R/bioconductor package to infer the hierarchical structure of cell populations from single-cell RNA-seq data. BMC Bioinformatics 17.](https://doi.org/10.1186/s12859-016-1175-6)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_celltree_maptpx <- create_ti_method_chooser(ti_celltree_maptpx, 'dynverse/celltree_maptpx')




#' Inferring a trajectory inference using [cellTree with vem](https://doi.org/10.1186/s12859-016-1175-6)
#' 
#' Will generate a trajectory using [cellTree with vem](https://doi.org/10.1186/s12859-016-1175-6). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/celltree_vem).
#' 
#' This methods was first wrapped inside R, see [ti_celltree_vem]
#' 
#' 
#' 
#' The method is described in: [duVerle, D.A., Yotsukura, S., Nomura, S., Aburatani, H., Tsuda, K., 2016. CellTree: an R/bioconductor package to infer the hierarchical structure of cell populations from single-cell RNA-seq data. BMC Bioinformatics 17.](https://doi.org/10.1186/s12859-016-1175-6)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_celltree_vem <- create_ti_method_chooser(ti_celltree_vem, 'dynverse/celltree_vem')




#' Inferring a trajectory inference using [DPT](https://doi.org/10.1038/nmeth.3971)
#' 
#' Will generate a trajectory using [DPT](https://doi.org/10.1038/nmeth.3971). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/dpt).
#' 
#' This methods was first wrapped inside R, see [ti_dpt]
#' 
#' The original code of this method is available [here](https://bioconductor.org/packages/release/bioc/html/destiny.html).
#' 
#' The method is described in: [Haghverdi, L., Büttner, M., Wolf, F.A., Buettner, F., Theis, F.J., 2016. Diffusion pseudotime robustly reconstructs lineage branching. Nature Methods 13, 845–848.](https://doi.org/10.1038/nmeth.3971)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_dpt <- create_ti_method_chooser(ti_dpt, 'dynverse/dpt')




#' Inferring a trajectory inference using [ElPiGraph cyclic](https://doi.org/https://github.com/Albluca/ElPiGraph.R)
#' 
#' Will generate a trajectory using [ElPiGraph cyclic](https://doi.org/https://github.com/Albluca/ElPiGraph.R). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/elpicycle).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/Albluca/ElPiGraph.R).
#' 
#' 
#' 
#' @param topology The kind of topology to detect \cr 
#'     discrete; default: "cycle"; possible values: cycle
#' @param NumNodes The number of nodes of the principal graph \cr 
#'     integer; default: 50L; possible values between 2 and 1000
#' @param NumEdges The maximum number of edges \cr 
#'     integer; default: 100000L; possible values between 2 and 100000
#' @param InitNodes Number of points to include in the initial graph \cr 
#'     integer; default: 2L; possible values between 2 and 1000
#' @param Mu Controls the elastic energy \cr 
#'     numeric; default: 0.1; possible values between 0.001 and 1
#' @param Lambda Controls the elastic energy \cr 
#'     numeric; default: 0.01; possible values between 0.001 and 1
#' @param MaxNumberOfIterations Maximum number of steps to embed the nodes \cr 
#'     integer; default: 10L; possible values between 1 and 1000
#' @param eps Minimal relative change in the position of the nodes to stop embedment \cr 
#'     numeric; default: 0.01; possible values between 0.001 and 1
#' @param CenterData Should data and initial node positions be centered? \cr 
#' 
#' @return The trajectory model
#' @export
ti_elpicycle <- function(
    topology = "cycle",
    NumNodes = 50L,
    NumEdges = 100000L,
    InitNodes = 2L,
    Mu = 0.1,
    Lambda = 0.01,
    MaxNumberOfIterations = 10L,
    eps = 0.01,
    CenterData = FALSE
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/elpicycle')
  do.call(method, args)
}




#' Inferring a trajectory inference using [ElPiGraph](https://doi.org/https://github.com/Albluca/ElPiGraph.R)
#' 
#' Will generate a trajectory using [ElPiGraph](https://doi.org/https://github.com/Albluca/ElPiGraph.R). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/elpigraph).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/Albluca/ElPiGraph.R).
#' 
#' 
#' 
#' @param topology The kind of topology to detect \cr 
#'     discrete; default: "tree"; possible values: cycle, tree, linear
#' @param NumNodes The number of nodes of the principal graph \cr 
#'     integer; default: 50L; possible values between 2 and 1000
#' @param NumEdges The maximum number of edges \cr 
#'     integer; default: 100000L; possible values between 2 and 100000
#' @param InitNodes Number of points to include in the initial graph \cr 
#'     integer; default: 2L; possible values between 2 and 1000
#' @param Mu Controls the elastic energy \cr 
#'     numeric; default: 0.1; possible values between 0.001 and 1
#' @param Lambda Controls the elastic energy \cr 
#'     numeric; default: 0.01; possible values between 0.001 and 1
#' @param MaxNumberOfIterations Maximum number of steps to embed the nodes \cr 
#'     integer; default: 10L; possible values between 1 and 1000
#' @param eps Minimal relative change in the position of the nodes to stop embedment \cr 
#'     numeric; default: 0.01; possible values between 0.001 and 1
#' @param CenterData Should data and initial node positions be centered? \cr 
#' 
#' @return The trajectory model
#' @export
ti_elpigraph <- function(
    topology = "tree",
    NumNodes = 50L,
    NumEdges = 100000L,
    InitNodes = 2L,
    Mu = 0.1,
    Lambda = 0.01,
    MaxNumberOfIterations = 10L,
    eps = 0.01,
    CenterData = FALSE
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/elpigraph')
  do.call(method, args)
}




#' Inferring a trajectory inference using [ElPiGraph linear](https://doi.org/https://github.com/Albluca/ElPiGraph.R)
#' 
#' Will generate a trajectory using [ElPiGraph linear](https://doi.org/https://github.com/Albluca/ElPiGraph.R). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/elpilinear).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/Albluca/ElPiGraph.R).
#' 
#' 
#' 
#' @param topology The kind of topology to detect \cr 
#'     discrete; default: "linear"; possible values: linear
#' @param NumNodes The number of nodes of the principal graph \cr 
#'     integer; default: 50L; possible values between 2 and 1000
#' @param NumEdges The maximum number of edges \cr 
#'     integer; default: 100000L; possible values between 2 and 100000
#' @param InitNodes Number of points to include in the initial graph \cr 
#'     integer; default: 2L; possible values between 2 and 1000
#' @param Mu Controls the elastic energy \cr 
#'     numeric; default: 0.1; possible values between 0.001 and 1
#' @param Lambda Controls the elastic energy \cr 
#'     numeric; default: 0.01; possible values between 0.001 and 1
#' @param MaxNumberOfIterations Maximum number of steps to embed the nodes \cr 
#'     integer; default: 10L; possible values between 1 and 1000
#' @param eps Minimal relative change in the position of the nodes to stop embedment \cr 
#'     numeric; default: 0.01; possible values between 0.001 and 1
#' @param CenterData Should data and initial node positions be centered? \cr 
#' 
#' @return The trajectory model
#' @export
ti_elpilinear <- function(
    topology = "linear",
    NumNodes = 50L,
    NumEdges = 100000L,
    InitNodes = 2L,
    Mu = 0.1,
    Lambda = 0.01,
    MaxNumberOfIterations = 10L,
    eps = 0.01,
    CenterData = FALSE
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/elpilinear')
  do.call(method, args)
}




#' Inferring a trajectory inference using [Embeddr](https://doi.org/10.1101/027219)
#' 
#' Will generate a trajectory using [Embeddr](https://doi.org/10.1101/027219). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/embeddr).
#' 
#' This methods was first wrapped inside R, see [ti_embeddr]
#' 
#' The original code of this method is available [here](https://github.com/kieranrcampbell/embeddr).
#' 
#' The method is described in: [Campbell, K., Ponting, C.P., Webber, C., 2015. Laplacian eigenmaps and principal curves for high resolution pseudotemporal ordering of single-cell RNA-seq profiles.](https://doi.org/10.1101/027219)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_embeddr <- create_ti_method_chooser(ti_embeddr, 'dynverse/embeddr')




#' Inferring a trajectory inference using Control: error
#' 
#' Will generate a trajectory using Control: error. This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/error).
#' 
#' This methods was first wrapped inside R, see [ti_error]
#' 
#' 
#' 
#' 
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_error <- create_ti_method_chooser(ti_error, 'dynverse/error')




#' Inferring a trajectory inference using [FateID](https://doi.org/10.1038/nmeth.4662)
#' 
#' Will generate a trajectory using [FateID](https://doi.org/10.1038/nmeth.4662). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/fateid).
#' 
#' 
#' 
#' The original code of this method is available [here](https://git.embl.de/velten/STEMNET).
#' 
#' The method is described in: [Herman, J.S., Sagar, Grün, D., 2018. FateID infers cell fate bias in multipotent progenitors from single-cell RNA-seq data. Nature Methods 15, 379–386.](https://doi.org/10.1038/nmeth.4662)
#' 
#' @param reclassify Whether to reclassify the cell grouping \cr 
#' @param clthr Real number between zero and one. This is the threshold for the fraction of random forest votes required to assign a cell not contained within the target clusters to one of these clusters. The value of this parameter should be sufficiently high to only reclassify cells with a high-confidence assignment. Default value is 0.9. \cr 
#'     numeric; default: 0.9; possible values between 0.1 and 1
#' @param nbfactor Positive integer number. Determines the number of trees grown for each random forest. The number of trees is given by the number of columns of th training set multiplied by \code{nbfactor}. Default value is 5. \cr 
#'     integer; default: 5L; possible values between 2 and 100
#' @param q Q real value between zero and one. This number specifies a threshold used for feature selection based on importance sampling. A reduced expression table is generated containing only features with an importance larger than the q-quantile for at least one of the classes (i. e. target clusters). Default value is 0.75. \cr 
#'     numeric; default: 0.75; possible values between 0 and 1
#' @param k Number of dimensions \cr 
#'     integer; default: 3L; possible values between 2 and 100
#' @param m Dimensionality reduction method to use. Can be tsne, cmd, dm or lle \cr 
#'     discrete; default: "tsne"; possible values: tsne, cmd, dm, lle
#' @param minnr Integer number of cells per target cluster to be selected for classification (test set) in each round of training. For each target cluster, the \code{minnr} cells with the highest similarity to a cell in the training set are selected for classification. If \code{z} is not \code{NULL} it is used as the similarity matrix for this step. Otherwise, \code{1-cor(x)} is used. Default value is 5. \cr 
#'     integer; default: 5L; possible values between 2 and 100
#' @param minnrh Integer number of cells from the training set used for classification. From each training set, the \code{minnrh} cells with the highest similarity to the training set are selected. If \code{z} is not \code{NULL} it is used as the similarity matrix for this step. Default value is 10. \cr 
#'     integer; default: 10L; possible values between 2 and 100
#' @param trthr Real value representing the threshold of the fraction of random forest votes required for the inclusion of a given cell for the computation of the principal curve. If \code{NULL} then only cells with a significant bias >1 are included for each trajectory. The bias is computed as the ratio of the number of votes for a trajectory and the number of votes for the trajectory with the second largest number of votes. By this means only the trajectory with the largest number of votes will receive a bias >1. The siginifcance is computed based on counting statistics on the difference in the number of votes. A significant bias requires a p-value < 0.05. \cr 
#'     numeric; default: 0.4; possible values between 0 and 1
#' 
#' @return The trajectory model
#' @export
ti_fateid <- function(
    reclassify = TRUE,
    clthr = 0.9,
    nbfactor = 5L,
    q = 0.75,
    k = 3L,
    m = "tsne",
    minnr = 5L,
    minnrh = 10L,
    trthr = 0.4
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/fateid')
  do.call(method, args)
}




#' Inferring a trajectory inference using [Growing Neural Gas](https://doi.org/https://github.com/rcannood/gng)
#' 
#' Will generate a trajectory using [Growing Neural Gas](https://doi.org/https://github.com/rcannood/gng). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/gng).
#' 
#' This methods was first wrapped inside R, see [ti_gng]
#' 
#' The original code of this method is available [here](https://github.com/rcannood/gng).
#' 
#' 
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_gng <- create_ti_method_chooser(ti_gng, 'dynverse/gng')




#' Inferring a trajectory inference using [GPfates](https://doi.org/10.1126/sciimmunol.aal2192)
#' 
#' Will generate a trajectory using [GPfates](https://doi.org/10.1126/sciimmunol.aal2192). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/gpfates).
#' 
#' 
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




#' Inferring a trajectory inference using [GrandPrix](https://doi.org/10.1101/227843)
#' 
#' Will generate a trajectory using [GrandPrix](https://doi.org/10.1101/227843). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/grandprix).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/ManchesterBioinference/GrandPrix).
#' 
#' The method is described in: [Ahmed, S., Rattray, M., Boukouvalas, A., 2017. GrandPrix: Scaling up the Bayesian GPLVM for single-cell data.](https://doi.org/10.1101/227843)
#' 
#' @param n_latent_dims  \cr 
#'     integer; default: 2L; possible values between 1 and 10
#' @param n_inducing_points  \cr 
#'     integer; default: 40L; possible values between 10 and 500
#' @param latent_prior_var  \cr 
#' @param latent_var  \cr 
#' 
#' @return The trajectory model
#' @export
ti_grandprix <- function(
    n_latent_dims = 2L,
    n_inducing_points = 40L,
    latent_prior_var = 0.1,
    latent_var = 0.028
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/grandprix')
  do.call(method, args)
}




#' Inferring a trajectory inference using Control: identity
#' 
#' Will generate a trajectory using Control: identity. This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/identity).
#' 
#' This methods was first wrapped inside R, see [ti_identity]
#' 
#' 
#' 
#' 
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_identity <- create_ti_method_chooser(ti_identity, 'dynverse/identity')




#' Inferring a trajectory inference using [MATCHER](https://doi.org/10.1186/s13059-017-1269-0)
#' 
#' Will generate a trajectory using [MATCHER](https://doi.org/10.1186/s13059-017-1269-0). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/matcher).
#' 
#' 
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




#' Inferring a trajectory inference using [MERLoT](https://doi.org/10.1101/261768)
#' 
#' Will generate a trajectory using [MERLoT](https://doi.org/10.1101/261768). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/merlot).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/soedinglab/merlot).
#' 
#' The method is described in: [Parra, R.G., Papadopoulos, N., Ahumada-Arranz, L., El Kholtei, J., Treutlein, B., Soeding, J., 2018. Reconstructing complex lineage trees from scRNA-seq data using MERLoT.](https://doi.org/10.1101/261768)
#' 
#' @param sigma Diffusion scale parameter of the Gaussian kernel. A larger sigma might be necessary if the eigenvalues can not be found because of a singularity in the matrix. Must be one of:\itemize{\item A character vector: \code{"local"} (default) or \code{"global"},\item a numeric global sigma -- a global sigma will be calculated using \code{\link[destiny:find_sigmas]{destiny::find_sigmas()}}\item or a \code{\link[destiny:Sigmas-class]{destiny::Sigmas-class()}} object.} \cr 
#'     discrete; default: "local"; possible values: local, global
#' @param distance A \code{\link[stats:dist]{stats::dist()}} object, or a character vector specifying which distance metric to use. Allowed measures:\itemize{\item Euclidean distance (default),\item cosine distance (1-corr(c_1, c_2)), or\item rank correlation distance (1-corr(rank(c_1), rank(c_2)))} \cr 
#'     discrete; default: "euclidean"; possible values: euclidean, cosine, rankcor
#' @param ndim Number of eigenvectors/dimensions to return \cr 
#'     integer; default: 20L; possible values between 2 and 20
#' @param density_norm Logical. If TRUE, use density normalisation \cr 
#'     logical; default: TRUE; possible values: TRUE, FALSE
#' @param n_local_lower If sigma == 'local', the \code{n_local_lower}:\code{n_local_upper} nearest neighbor(s) determine(s) the local sigma \cr 
#'     integer; default: 5L; possible values between 2 and 20
#' @param n_local_upper See \code{n_local_lower} \cr 
#'     integer; default: 7L; possible values between 2 and 20
#' @param w_width Window width to use for deciding the branch cutoff \cr 
#'     numeric; default: 0.01; possible values between 1e-04 and 1
#' @param n_components_to_use Which components to use in downstream analysis \cr 
#'     integer; default: 3; possible values between 2 and 20
#' @param N_yk Number of nodes for the elastic principal tree \cr 
#'     integer; default: 100; possible values between 2 and 1000
#' @param lambda_0 Principal elastic tree energy function parameter. \cr 
#'     numeric; default: 1e-10; possible values between 1e-15 and 1e-04
#' @param mu_0 Principal elastic tree energy function parameter. \cr 
#'     numeric; default: 0.0025; possible values between 5e-04 and 0.005
#' @param increaseFactor_mu Factor by which the mu will be increased for the embedding \cr 
#'     numeric; default: 20; possible values between 2 and 50
#' @param increaseFactor_lambda Factor by which the mu will be increased for the embedding \cr 
#'     numeric; default: 20; possible values between 2 and 50
#' @param FixEndpoints Documentation not provided by authors \cr 
#'     logical; default: FALSE; possible values: TRUE, FALSE
#' 
#' @return The trajectory model
#' @export
ti_merlot <- function(
    sigma = "local",
    distance = "euclidean",
    ndim = 20L,
    density_norm = TRUE,
    n_local_lower = 5L,
    n_local_upper = 7L,
    w_width = 0.01,
    n_components_to_use = 3,
    N_yk = 100,
    lambda_0 = 1e-10,
    mu_0 = 0.0025,
    increaseFactor_mu = 20,
    increaseFactor_lambda = 20,
    FixEndpoints = FALSE
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/merlot')
  do.call(method, args)
}




#' Inferring a trajectory inference using [mfa](https://doi.org/10.12688/wellcomeopenres.11087.1)
#' 
#' Will generate a trajectory using [mfa](https://doi.org/10.12688/wellcomeopenres.11087.1). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/mfa).
#' 
#' This methods was first wrapped inside R, see [ti_mfa]
#' 
#' The original code of this method is available [here](https://github.com/kieranrcampbell/mfa).
#' 
#' The method is described in: [Campbell, K.R., Yau, C., 2017. Probabilistic modeling of bifurcations in single-cell gene expression data using a Bayesian mixture of factor analyzers. Wellcome Open Research 2, 19.](https://doi.org/10.12688/wellcomeopenres.11087.1)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_mfa <- create_ti_method_chooser(ti_mfa, 'dynverse/mfa')




#' Inferring a trajectory inference using [Monocle DDRTree](https://doi.org/10.1038/nmeth.4402)
#' 
#' Will generate a trajectory using [Monocle DDRTree](https://doi.org/10.1038/nmeth.4402). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/monocle_ddrtree).
#' 
#' This methods was first wrapped inside R, see [ti_monocle_ddrtree]
#' 
#' The original code of this method is available [here](https://github.com/cole-trapnell-lab/monocle-release).
#' 
#' The method is described in: [Qiu, X., Mao, Q., Tang, Y., Wang, L., Chawla, R., Pliner, H.A., Trapnell, C., 2017. Reversed graph embedding resolves complex single-cell trajectories. Nature Methods 14, 979–982.](https://doi.org/10.1038/nmeth.4402)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_monocle_ddrtree <- create_ti_method_chooser(ti_monocle_ddrtree, 'dynverse/monocle_ddrtree')




#' Inferring a trajectory inference using [Monocle ICA](https://doi.org/10.1038/nmeth.4402)
#' 
#' Will generate a trajectory using [Monocle ICA](https://doi.org/10.1038/nmeth.4402). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/monocle_ica).
#' 
#' This methods was first wrapped inside R, see [ti_monocle_ica]
#' 
#' The original code of this method is available [here](https://github.com/cole-trapnell-lab/monocle-release).
#' 
#' The method is described in: [Qiu, X., Mao, Q., Tang, Y., Wang, L., Chawla, R., Pliner, H.A., Trapnell, C., 2017. Reversed graph embedding resolves complex single-cell trajectories. Nature Methods 14, 979–982.](https://doi.org/10.1038/nmeth.4402)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_monocle_ica <- create_ti_method_chooser(ti_monocle_ica, 'dynverse/monocle_ica')




#' Inferring a trajectory inference using [Mpath](https://doi.org/10.1038/ncomms11988)
#' 
#' Will generate a trajectory using [Mpath](https://doi.org/10.1038/ncomms11988). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/mpath).
#' 
#' This methods was first wrapped inside R, see [ti_mpath]
#' 
#' The original code of this method is available [here](https://github.com/JinmiaoChenLab/Mpath).
#' 
#' The method is described in: [Chen, J., Schlitzer, A., Chakarov, S., Ginhoux, F., Poidinger, M., 2016. Mpath maps multi-branching single-cell trajectories revealing progenitor cell progression during development. Nature Communications 7, 11988.](https://doi.org/10.1038/ncomms11988)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_mpath <- create_ti_method_chooser(ti_mpath, 'dynverse/mpath')




#' Inferring a trajectory inference using [ouija](https://doi.org/10.1101/060442)
#' 
#' Will generate a trajectory using [ouija](https://doi.org/10.1101/060442). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/ouija).
#' 
#' This methods was first wrapped inside R, see [ti_ouija]
#' 
#' The original code of this method is available [here](https://github.com/kieranrcampbell/ouija).
#' 
#' The method is described in: [Campbell, K.R., Yau, C., 2016. A descriptive marker gene approach to single-cell pseudotime inference.](https://doi.org/10.1101/060442)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_ouija <- create_ti_method_chooser(ti_ouija, 'dynverse/ouija')




#' Inferring a trajectory inference using [ouijaflow](https://doi.org/10.1101/060442)
#' 
#' Will generate a trajectory using [ouijaflow](https://doi.org/10.1101/060442). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/ouijaflow).
#' 
#' 
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




#' Inferring a trajectory inference using [PAGA](https://doi.org/10.1101/208819)
#' 
#' Will generate a trajectory using [PAGA](https://doi.org/10.1101/208819). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/paga).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/theislab/graph_abstraction).
#' 
#' The method is described in: [Wolf, F.A., Hamey, F., Plass, M., Solana, J., Dahlin, J.S., Gottgens, B., Rajewsky, N., Simon, L., Theis, F.J., 2017. Graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells.](https://doi.org/10.1101/208819)
#' 
#' @param n_neighbors Number of neighbours for knn \cr 
#'     integer; default: 30L; possible values between 1 and 100
#' @param n_comps Number of principal components \cr 
#'     integer; default: 50L; possible values between 0 and 100
#' @param n_dcs Number of diffusion components for denoising graph, 0 means no denoising. \cr 
#'     integer; default: 15L; possible values between 0 and 40
#' @param resolution Resolution of louvain clustering, which determines the granularity of the clustering. Higher values will result in more clusters. \cr 
#'     numeric; default: 2.5; possible values between 0.1 and 10
#' @param embedding_type Either 'umap' (scales very well, recommended for very large datasets) or 'fa' (ForceAtlas2, often a bit more intuitive for small datasets). \cr 
#' 
#' @return The trajectory model
#' @export
ti_paga <- function(
    n_neighbors = 30L,
    n_comps = 50L,
    n_dcs = 15L,
    resolution = 2.5,
    embedding_type = "fa"
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/paga')
  do.call(method, args)
}




#' Inferring a trajectory inference using [pCreode](https://doi.org/10.1016/j.cels.2017.10.012)
#' 
#' Will generate a trajectory using [pCreode](https://doi.org/10.1016/j.cels.2017.10.012). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/pcreode).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/KenLauLab/pCreode).
#' 
#' The method is described in: [Herring, C.A., Banerjee, A., McKinley, E.T., Simmons, A.J., Ping, J., Roland, J.T., Franklin, J.L., Liu, Q., Gerdes, M.J., Coffey, R.J., Lau, K.S., 2018. Unsupervised Trajectory Analysis of Single-Cell RNA-Seq and Imaging Data Reveals Alternative Tuft Cell Origins in the Gut. Cell Systems 6, 37–51.e9.](https://doi.org/10.1016/j.cels.2017.10.012)
#' 
#' @param n_pca_components  \cr 
#'     integer; default: 3L; possible values between 2 and 10
#' @param radius  \cr 
#'     numeric; default: 1L; possible values between 0.01 and 10
#' @param noise  \cr 
#'     numeric; default: 8L; possible values between 1 and 20
#' @param target  \cr 
#'     numeric; default: 25L; possible values between 5 and 100
#' @param num_runs  \cr 
#'     integer; default: 10L; possible values between 10 and 1000
#' 
#' @return The trajectory model
#' @export
ti_pcreode <- function(
    n_pca_components = 3L,
    radius = 1L,
    noise = 8L,
    target = 25L,
    num_runs = 10L
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/pcreode')
  do.call(method, args)
}




#' Inferring a trajectory inference using Periodic PrinCurve
#' 
#' Will generate a trajectory using Periodic PrinCurve. This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/periodpc).
#' 
#' This methods was first wrapped inside R, see [ti_periodpc]
#' 
#' 
#' 
#' 
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_periodpc <- create_ti_method_chooser(ti_periodpc, 'dynverse/periodpc')




#' Inferring a trajectory inference using [PhenoPath](https://doi.org/10.1101/159913)
#' 
#' Will generate a trajectory using [PhenoPath](https://doi.org/10.1101/159913). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/phenopath).
#' 
#' This methods was first wrapped inside R, see [ti_phenopath]
#' 
#' The original code of this method is available [here](https://github.com/kieranrcampbell/phenopath).
#' 
#' The method is described in: [Campbell, K., Yau, C., 2017. Uncovering genomic trajectories with heterogeneous genetic and environmental backgrounds across single-cells and populations.](https://doi.org/10.1101/159913)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_phenopath <- create_ti_method_chooser(ti_phenopath, 'dynverse/phenopath')




#' Inferring a trajectory inference using [projected PAGA](https://doi.org/10.1101/208819)
#' 
#' Will generate a trajectory using [projected PAGA](https://doi.org/10.1101/208819). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/praga).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/theislab/graph_abstraction).
#' 
#' The method is described in: [Wolf, F.A., Hamey, F., Plass, M., Solana, J., Dahlin, J.S., Gottgens, B., Rajewsky, N., Simon, L., Theis, F.J., 2017. Graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells.](https://doi.org/10.1101/208819)
#' 
#' @param embedding_type Either 'umap' (scales very well, recommended for very large datasets) or 'fa' (ForceAtlas2, often a bit more intuitive for small datasets). \cr 
#' @param n_comps Number of principal components \cr 
#'     integer; default: 50L; possible values between 0 and 100
#' @param n_dcs Number of diffusion components for denoising graph, 0 means no denoising. \cr 
#'     integer; default: 15L; possible values between 0 and 40
#' @param n_neighbors Number of neighbours for knn \cr 
#'     integer; default: 30L; possible values between 1 and 100
#' @param resolution Resolution of louvain clustering, which determines the granularity of the clustering. Higher values will result in more clusters. \cr 
#'     numeric; default: 2.5; possible values between 0.1 and 10
#' 
#' @return The trajectory model
#' @export
ti_praga <- function(
    embedding_type = "fa",
    n_comps = 50L,
    n_dcs = 15L,
    n_neighbors = 30L,
    resolution = 2.5
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/praga')
  do.call(method, args)
}




#' Inferring a trajectory inference using [pseudogp](https://doi.org/10.1371/journal.pcbi.1005212)
#' 
#' Will generate a trajectory using [pseudogp](https://doi.org/10.1371/journal.pcbi.1005212). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/pseudogp).
#' 
#' This methods was first wrapped inside R, see [ti_pseudogp]
#' 
#' The original code of this method is available [here](https://github.com/kieranrcampbell/pseudogp).
#' 
#' The method is described in: [Campbell, K.R., Yau, C., 2016. Order Under Uncertainty: Robust Differential Expression Analysis Using Probabilistic Models for Pseudotime Inference. PLOS Computational Biology 12, e1005212.](https://doi.org/10.1371/journal.pcbi.1005212)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_pseudogp <- create_ti_method_chooser(ti_pseudogp, 'dynverse/pseudogp')




#' Inferring a trajectory inference using Control: random
#' 
#' Will generate a trajectory using Control: random. This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/random).
#' 
#' This methods was first wrapped inside R, see [ti_random]
#' 
#' 
#' 
#' 
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_random <- create_ti_method_chooser(ti_random, 'dynverse/random')




#' Inferring a trajectory inference using [reCAT](https://doi.org/10.1038/s41467-017-00039-z)
#' 
#' Will generate a trajectory using [reCAT](https://doi.org/10.1038/s41467-017-00039-z). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/recat).
#' 
#' This methods was first wrapped inside R, see [ti_recat]
#' 
#' The original code of this method is available [here](https://github.com/tinglab/reCAT).
#' 
#' The method is described in: [Liu, Z., Lou, H., Xie, K., Wang, H., Chen, N., Aparicio, O.M., Zhang, M.Q., Jiang, R., Chen, T., 2017. Reconstructing cell cycle pseudo time-series via single-cell transcriptome data. Nature Communications 8.](https://doi.org/10.1038/s41467-017-00039-z)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_recat <- create_ti_method_chooser(ti_recat, 'dynverse/recat')




#' Inferring a trajectory inference using [SCIMITAR](https://doi.org/10.1142/9789813207813_0053)
#' 
#' Will generate a trajectory using [SCIMITAR](https://doi.org/10.1142/9789813207813_0053). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/scimitar).
#' 
#' 
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




#' Inferring a trajectory inference using [SCORPIUS](https://doi.org/10.1101/079509)
#' 
#' Will generate a trajectory using [SCORPIUS](https://doi.org/10.1101/079509). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/scorpius).
#' 
#' This methods was first wrapped inside R, see [ti_scorpius]
#' 
#' The original code of this method is available [here](https://github.com/rcannood/SCORPIUS).
#' 
#' The method is described in: [Cannoodt, R., Saelens, W., Sichien, D., Tavernier, S., Janssens, S., Guilliams, M., Lambrecht, B.N., De Preter, K., Saeys, Y., 2016. SCORPIUS improves trajectory inference and identifies novel modules in dendritic cell development.](https://doi.org/10.1101/079509)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_scorpius <- create_ti_method_chooser(ti_scorpius, 'dynverse/scorpius')




#' Inferring a trajectory inference using [SCORPIUS sparse](https://doi.org/10.1101/079509)
#' 
#' Will generate a trajectory using [SCORPIUS sparse](https://doi.org/10.1101/079509). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/scorpius_sparse).
#' 
#' This methods was first wrapped inside R, see [ti_scorpius_sparse]
#' 
#' The original code of this method is available [here](https://github.com/rcannood/SCORPIUS).
#' 
#' The method is described in: [Cannoodt, R., Saelens, W., Sichien, D., Tavernier, S., Janssens, S., Guilliams, M., Lambrecht, B.N., De Preter, K., Saeys, Y., 2016. SCORPIUS improves trajectory inference and identifies novel modules in dendritic cell development.](https://doi.org/10.1101/079509)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_scorpius_sparse <- create_ti_method_chooser(ti_scorpius_sparse, 'dynverse/scorpius_sparse')




#' Inferring a trajectory inference using [SCOUP](https://doi.org/10.1186/s12859-016-1109-3)
#' 
#' Will generate a trajectory using [SCOUP](https://doi.org/10.1186/s12859-016-1109-3). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/scoup).
#' 
#' This methods was first wrapped inside R, see [ti_scoup]
#' 
#' The original code of this method is available [here](https://github.com/gcyuan/PySCUBA).
#' 
#' The method is described in: [Matsumoto, H., Kiryu, H., 2016. SCOUP: a probabilistic model based on the Ornstein–Uhlenbeck process to analyze single-cell expression data during differentiation. BMC Bioinformatics 17.](https://doi.org/10.1186/s12859-016-1109-3)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_scoup <- create_ti_method_chooser(ti_scoup, 'dynverse/scoup')




#' Inferring a trajectory inference using [SCUBA](https://doi.org/10.1073/pnas.1408993111)
#' 
#' Will generate a trajectory using [SCUBA](https://doi.org/10.1073/pnas.1408993111). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/scuba).
#' 
#' 
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




#' Inferring a trajectory inference using Control: shuffle
#' 
#' Will generate a trajectory using Control: shuffle. This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/shuffle).
#' 
#' This methods was first wrapped inside R, see [ti_shuffle]
#' 
#' 
#' 
#' 
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_shuffle <- create_ti_method_chooser(ti_shuffle, 'dynverse/shuffle')




#' Inferring a trajectory inference using [Sincell](https://doi.org/10.1093/bioinformatics/btv368)
#' 
#' Will generate a trajectory using [Sincell](https://doi.org/10.1093/bioinformatics/btv368). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/sincell).
#' 
#' This methods was first wrapped inside R, see [ti_sincell]
#' 
#' The original code of this method is available [here](https://github.com/Cortalak/MCA_Sincell_0).
#' 
#' The method is described in: [Juliá, M., Telenti, A., Rausell, A., 2015. Sincell: an R/Bioconductor package for statistical assessment of cell-state hierarchies from single-cell RNA-seq: Fig. 1. Bioinformatics 31, 3380–3382.](https://doi.org/10.1093/bioinformatics/btv368)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_sincell <- create_ti_method_chooser(ti_sincell, 'dynverse/sincell')




#' Inferring a trajectory inference using [SLICE](https://doi.org/10.1093/nar/gkw1278)
#' 
#' Will generate a trajectory using [SLICE](https://doi.org/10.1093/nar/gkw1278). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/slice).
#' 
#' This methods was first wrapped inside R, see [ti_slice]
#' 
#' The original code of this method is available [here](https://research.cchmc.org/pbge/slice.html).
#' 
#' The method is described in: [Guo, M., Bao, E.L., Wagner, M., Whitsett, J.A., Xu, Y., 2016. SLICE: determining cell differentiation and lineage based on single cell entropy. Nucleic Acids Research gkw1278.](https://doi.org/10.1093/nar/gkw1278)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_slice <- create_ti_method_chooser(ti_slice, 'dynverse/slice')




#' Inferring a trajectory inference using [SLICER](https://doi.org/10.1186/s13059-016-0975-3)
#' 
#' Will generate a trajectory using [SLICER](https://doi.org/10.1186/s13059-016-0975-3). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/slicer).
#' 
#' This methods was first wrapped inside R, see [ti_slicer]
#' 
#' The original code of this method is available [here](https://github.com/jw156605/SLICER).
#' 
#' The method is described in: [Welch, J.D., Hartemink, A.J., Prins, J.F., 2016. SLICER: inferring branched, nonlinear cellular trajectories from single cell RNA-seq data. Genome Biology 17.](https://doi.org/10.1186/s13059-016-0975-3)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_slicer <- create_ti_method_chooser(ti_slicer, 'dynverse/slicer')




#' Inferring a trajectory inference using [Slingshot](https://doi.org/10.1101/128843)
#' 
#' Will generate a trajectory using [Slingshot](https://doi.org/10.1101/128843). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/slingshot).
#' 
#' This methods was first wrapped inside R, see [ti_slingshot]
#' 
#' The original code of this method is available [here](https://github.com/kstreet13/slingshot).
#' 
#' The method is described in: [Street, K., Risso, D., Fletcher, R.B., Das, D., Ngai, J., Yosef, N., Purdom, E., Dudoit, S., 2017. Slingshot: Cell lineage and pseudotime inference for single-cell transcriptomics.](https://doi.org/10.1101/128843)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_slingshot <- create_ti_method_chooser(ti_slingshot, 'dynverse/slingshot')




#' Inferring a trajectory inference using [StemID](https://doi.org/10.1016/j.stem.2016.05.010)
#' 
#' Will generate a trajectory using [StemID](https://doi.org/10.1016/j.stem.2016.05.010). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/stemid).
#' 
#' This methods was first wrapped inside R, see [ti_stemid]
#' 
#' The original code of this method is available [here](https://github.com/dgrun/StemID).
#' 
#' The method is described in: [Grün, D., Muraro, M.J., Boisset, J.-C., Wiebrands, K., Lyubimova, A., Dharmadhikari, G., van den Born, M., van Es, J., Jansen, E., Clevers, H., de Koning, E.J.P., van Oudenaarden, A., 2016. De Novo Prediction of Stem Cell Identity using Single-Cell Transcriptome Data. Cell Stem Cell 19, 266–277.](https://doi.org/10.1016/j.stem.2016.05.010)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_stemid <- create_ti_method_chooser(ti_stemid, 'dynverse/stemid')




#' Inferring a trajectory inference using [StemID2](https://doi.org/10.1016/j.stem.2016.05.010)
#' 
#' Will generate a trajectory using [StemID2](https://doi.org/10.1016/j.stem.2016.05.010). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/stemid2).
#' 
#' This methods was first wrapped inside R, see [ti_stemid2]
#' 
#' The original code of this method is available [here](https://github.com/dgrun/RaceID3_StemID2).
#' 
#' The method is described in: [Grün, D., Muraro, M.J., Boisset, J.-C., Wiebrands, K., Lyubimova, A., Dharmadhikari, G., van den Born, M., van Es, J., Jansen, E., Clevers, H., de Koning, E.J.P., van Oudenaarden, A., 2016. De Novo Prediction of Stem Cell Identity using Single-Cell Transcriptome Data. Cell Stem Cell 19, 266–277.](https://doi.org/10.1016/j.stem.2016.05.010)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_stemid2 <- create_ti_method_chooser(ti_stemid2, 'dynverse/stemid2')




#' Inferring a trajectory inference using [STEMNET](https://doi.org/10.1038/ncb3493)
#' 
#' Will generate a trajectory using [STEMNET](https://doi.org/10.1038/ncb3493). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/stemnet).
#' 
#' 
#' 
#' The original code of this method is available [here](https://git.embl.de/velten/STEMNET).
#' 
#' The method is described in: [Velten, L., Haas, S.F., Raffel, S., Blaszkiewicz, S., Islam, S., Hennig, B.P., Hirche, C., Lutz, C., Buss, E.C., Nowak, D., Boch, T., Hofmann, W.-K., Ho, A.D., Huber, W., Trumpp, A., Essers, M.A.G., Steinmetz, L.M., 2017. Human haematopoietic stem cell lineage commitment is a continuous process. Nature Cell Biology 19, 271–281.](https://doi.org/10.1038/ncb3493)
#' 
#' @param alpha The elastic net mixing parameter of the ‘glmnet’ classifier. \cr 
#'     numeric; default: 0.1; possible values between 0.001 and 10
#' @param lambda_auto Whether to select the lambda by cross-validation \cr 
#' @param lambda The lambda penalty of GLM. \cr 
#'     numeric; default: 0.1; possible values between 0.05 and 1
#' 
#' @return The trajectory model
#' @export
ti_stemnet <- function(
    alpha = 0.1,
    lambda_auto = TRUE,
    lambda = 0.1
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/stemnet')
  do.call(method, args)
}




#' Inferring a trajectory inference using [topslam](https://doi.org/10.1101/057778)
#' 
#' Will generate a trajectory using [topslam](https://doi.org/10.1101/057778). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/topslam).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/mzwiessele/topslam).
#' 
#' The method is described in: [Zwiessele, M., Lawrence, N.D., 2016. Topslam: Waddington Landscape Recovery for Single Cell Experiments.](https://doi.org/10.1101/057778)
#' 
#' @param n_components The number of components \cr 
#'     integer; default: 2L; possible values between 2 and 10
#' @param n_neighbors The number of neighbors \cr 
#'     integer; default: 10L; possible values between 2 and 100
#' @param linear_dims  \cr 
#'     integer; default: 0L; possible values between 0 and 5
#' @param max_iters The number of iterations to optimize over \cr 
#'     integer; default: 1000L; possible values between 10 and 10000
#' @param dimreds Which dimensionality reductions to use; tSNE, PCA, Spectral, Isomap and/or ICA \cr 
#' 
#' @return The trajectory model
#' @export
ti_topslam <- function(
    n_components = 2L,
    n_neighbors = 10L,
    linear_dims = 0L,
    max_iters = 1000L,
    dimreds = c(TRUE, TRUE, TRUE, TRUE, TRUE)
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/topslam')
  do.call(method, args)
}




#' Inferring a trajectory inference using [TSCAN](https://doi.org/10.1093/nar/gkw430)
#' 
#' Will generate a trajectory using [TSCAN](https://doi.org/10.1093/nar/gkw430). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/tscan).
#' 
#' This methods was first wrapped inside R, see [ti_tscan]
#' 
#' The original code of this method is available [here](https://github.com/zji90/TSCAN).
#' 
#' The method is described in: [Ji, Z., Ji, H., 2016. TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis. Nucleic Acids Research 44, e117–e117.](https://doi.org/10.1093/nar/gkw430)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_tscan <- create_ti_method_chooser(ti_tscan, 'dynverse/tscan')




#' Inferring a trajectory inference using [wanderlust](https://doi.org/10.1038/nbt.3569)
#' 
#' Will generate a trajectory using [wanderlust](https://doi.org/10.1038/nbt.3569). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/wanderlust).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/ManuSetty/wishbone).
#' 
#' The method is described in: [Setty, M., Tadmor, M.D., Reich-Zeliger, S., Angel, O., Salame, T.M., Kathail, P., Choi, K., Bendall, S., Friedman, N., Pe’er, D., 2016. Wishbone identifies bifurcating developmental trajectories from single-cell data. Nature Biotechnology 34, 637–645.](https://doi.org/10.1038/nbt.3569)
#' 
#' @param branch Whether to allow a single bifurcation within the trajectory (wishbone versus wanderlust) \cr 
#' @param epsilon Epsilon \cr 
#'     numeric; default: 1L; possible values between 0.1 and 5
#' @param k K parameter \cr 
#'     integer; default: 25L; possible values between 15 and 100
#' @param knn K-nearest neighbours for diffusion \cr 
#'     integer; default: 25L; possible values between 15 and 100
#' @param n_diffusion_components Number of diffusion components \cr 
#'     integer; default: 3L; possible values between 3 and 20
#' @param n_pca_components Number of pca components \cr 
#'     integer; default: 30L; possible values between 15 and 100
#' @param normalise  \cr 
#' @param num_waypoints Number of waypoints \cr 
#'     integer; default: 250L; possible values between 100 and 500
#' 
#' @return The trajectory model
#' @export
ti_wanderlust <- function(
    branch = FALSE,
    epsilon = 1L,
    k = 25L,
    knn = 25L,
    n_diffusion_components = 3L,
    n_pca_components = 30L,
    normalise = TRUE,
    num_waypoints = 250L
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/wanderlust')
  do.call(method, args)
}




#' Inferring a trajectory inference using [Waterfall](https://doi.org/10.1016/j.stem.2015.07.013)
#' 
#' Will generate a trajectory using [Waterfall](https://doi.org/10.1016/j.stem.2015.07.013). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/waterfall).
#' 
#' This methods was first wrapped inside R, see [ti_waterfall]
#' 
#' The original code of this method is available [here](http://www.cell.com/cms/attachment/2038326541/2052521637/mmc9.zip).
#' 
#' The method is described in: [Shin, J., Berg, D.A., Zhu, Y., Shin, J.Y., Song, J., Bonaguidi, M.A., Enikolopov, G., Nauen, D.W., Christian, K.M., Ming, G., Song, H., 2015. Single-Cell RNA-Seq with Waterfall Reveals Molecular Cascades underlying Adult Neurogenesis. Cell Stem Cell 17, 360–372.](https://doi.org/10.1016/j.stem.2015.07.013)
#' 
#' @param docker Whether to use the docker container or the R wrapper
#' 
#' @return The trajectory model
#' @export
ti_waterfall <- create_ti_method_chooser(ti_waterfall, 'dynverse/waterfall')




#' Inferring a trajectory inference using [Wishbone](https://doi.org/10.1038/nbt.3569)
#' 
#' Will generate a trajectory using [Wishbone](https://doi.org/10.1038/nbt.3569). This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/wishbone).
#' 
#' 
#' 
#' The original code of this method is available [here](https://github.com/ManuSetty/wishbone).
#' 
#' The method is described in: [Setty, M., Tadmor, M.D., Reich-Zeliger, S., Angel, O., Salame, T.M., Kathail, P., Choi, K., Bendall, S., Friedman, N., Pe’er, D., 2016. Wishbone identifies bifurcating developmental trajectories from single-cell data. Nature Biotechnology 34, 637–645.](https://doi.org/10.1038/nbt.3569)
#' 
#' @param normalise  \cr 
#' @param knn K-nearest neighbours for diffusion \cr 
#'     integer; default: 25L; possible values between 15 and 100
#' @param n_diffusion_components Number of diffusion components \cr 
#'     integer; default: 3L; possible values between 3 and 20
#' @param n_pca_components Number of pca components \cr 
#'     integer; default: 30L; possible values between 15 and 100
#' @param k K parameter \cr 
#'     integer; default: 25L; possible values between 15 and 100
#' @param num_waypoints Number of waypoints \cr 
#'     integer; default: 250L; possible values between 100 and 500
#' @param epsilon Epsilon \cr 
#'     numeric; default: 1L; possible values between 0.1 and 5
#' @param branch Whether to allow a single bifurcation within the trajectory (wishbone versus wanderlust) \cr 
#' 
#' @return The trajectory model
#' @export
ti_wishbone <- function(
    normalise = TRUE,
    knn = 25L,
    n_diffusion_components = 3L,
    n_pca_components = 30L,
    k = 25L,
    num_waypoints = 250L,
    epsilon = 1L,
    branch = TRUE
) {
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/wishbone')
  do.call(method, args)
}



