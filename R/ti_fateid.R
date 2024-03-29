######################################### DO NOT EDIT! #########################################
#### This file is automatically generated from data-raw/2-generate_r_code_from_containers.R ####
################################################################################################

#' @title FateID
#' 
#' @description
#' Will generate a trajectory using [FateID](https://doi.org/10.1038/nmeth.4662).
#' 
#' This method was wrapped inside a
#' [container](https://github.com/dynverse/ti_fateid).
#' The original code of this method is available
#' [here](https://github.com/dgrun/FateID).
#' 
#' @references Herman, J.S., Sagar, Grün, D., 2018. FateID infers cell fate bias
#' in multipotent progenitors from single-cell RNA-seq data. Nature Methods 15,
#' 379–386.
#' 
#' @param reclassify Whether to reclassify the cell grouping. Default: TRUE.
#' Format: logical.
#' @param clthr Real number between zero and one. This is the threshold for the
#' fraction of random forest votes required to assign a cell not contained within
#' the target clusters to one of these clusters. The value of this parameter
#' should be sufficiently high to only reclassify cells with a high-confidence
#' assignment. Default value is 0.9. Domain: U(0.1, 1). Default: 0.9. Format:
#' numeric.
#' @param nbfactor Positive integer number. Determines the number of trees grown
#' for each random forest. The number of trees is given by the number of columns
#' of th training set multiplied by `nbfactor`. Default value is 5. Domain: U(2,
#' 100). Default: 5. Format: integer.
#' @param q Q real value between zero and one. This number specifies a threshold
#' used for feature selection based on importance sampling. A reduced expression
#' table is generated containing only features with an importance larger than the
#' q-quantile for at least one of the classes (i. e. target clusters). Default
#' value is 0.75. Domain: U(0, 1). Default: 0.75. Format: numeric.
#' @param k Number of dimensions. Domain: U(2, 100). Default: 3. Format: integer.
#' @param m Dimensionality reduction method to use. Can be tsne, cmd, dm or lle.
#' Domain: {tsne, cmd, dm, lle}. Default: tsne. Format: character.
#' @param minnr Integer number of cells per target cluster to be selected for
#' classification (test set) in each round of training. For each target cluster,
#' the `minnr` cells with the highest similarity to a cell in the training set are
#' selected for classification. If `z` is not `NULL` it is used as the similarity
#' matrix for this step. Otherwise, `1-cor(x)` is used. Default value is 5.
#' Domain: U(2, 100). Default: 5. Format: integer.
#' @param minnrh Integer number of cells from the training set used for
#' classification. From each training set, the `minnrh` cells with the highest
#' similarity to the training set are selected. If `z` is not `NULL` it is used as
#' the similarity matrix for this step. Default value is 10. Domain: U(2, 100).
#' Default: 10. Format: integer.
#' @param trthr Real value representing the threshold of the fraction of random
#' forest votes required for the inclusion of a given cell for the computation of
#' the principal curve. If `NULL` then only cells with a significant bias >1 are
#' included for each trajectory. The bias is computed as the ratio of the number
#' of votes for a trajectory and the number of votes for the trajectory with the
#' second largest number of votes. By this means only the trajectory with the
#' largest number of votes will receive a bias >1. The siginifcance is computed
#' based on counting statistics on the difference in the number of votes. A
#' significant bias requires a p-value < 0.05. Domain: U(0, 1). Default: 0.4.
#' Format: numeric.
#' @param force Do not use! This is a parameter to force FateID to run on
#' benchmark datasets where not enough end groups are present. Default: FALSE.
#' Format: logical.
#' 
#' @keywords method
#' 
#' @return A TI method wrapper to be used together with
#' \code{\link[dynwrap:infer_trajectories]{infer_trajectory}}
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
    trthr = 0.4,
    force = FALSE
) {
  method_choose_backend(
    package_repository = NULL,
    package_name = NULL,
    function_name = NULL,
    package_version = NULL,
    container_id = "dynverse/ti_fateid:v0.9.9.01"
  )(
    reclassify = reclassify,
    clthr = clthr,
    nbfactor = nbfactor,
    q = q,
    k = k,
    m = m,
    minnr = minnr,
    minnrh = minnrh,
    trthr = trthr,
    force = force
  )
}

