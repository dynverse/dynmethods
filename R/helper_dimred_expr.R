#' Perform simple dimensionality reduction
#'
#' @param x Log transformed expression data, with rows as cells and columns as features
#' @param method The name of the dimensionality reduction method to use
#' @param ndim The number of dimensions
#' @param ... Any arguments to be passed to the dimensionality reduction method
#'
#' @rdname dimred
#' @export
list_dimred_methods <- function() {
  list(
    pca = dimred_pca,
    mds = dimred_mds,
    tsne = dimred_tsne,
    ica = dimred_ica,
    lle = dimred_lle
  )
}

#' @rdname dimred
#' @export
dimred <- function(x, method, ndim, ...) {
  methods <- list_dimred_methods()
  if (method %in% names(methods)) {
    meth <- methods[[method]]
    params <- list(x = x, ndim = ndim, ...)
    do.call(meth, params)
  } else {
    stop("Method ", sQuote(method), " not found.")
  }
}

process_dimred <- function(space, rn) {
  space <- as.matrix(space)
  dimnames(space) <- list(rn, paste0("Comp", seq_len(ncol(space))))
  space
}

dimred_pca <- function(x, ndim = 3) {
  requireNamespace("stats")
  space <- stats::prcomp(x)$x[,seq_len(ndim)]
  process_dimred(space, rownames(x))
}

dimred_mds <- function(x, ndim = 3) {
  requireNamespace("SCORPIUS")
  space <- SCORPIUS::reduce_dimensionality(SCORPIUS::correlation_distance(x), ndim = ndim)
  process_dimred(space, rownames(x))
}

dimred_tsne <- function(x, ndim = 3) {
  requireNamespace("SCORPIUS")
  requireNamespace("Rtsne")
  requireNamespace("stats")
  space <- Rtsne::Rtsne(stats::as.dist(SCORPIUS::correlation_distance(x)), dims = ndim, is_distance = TRUE)$Y
  process_dimred(space, rownames(x))
}

dimred_ica <- function(x, ndim = 3) {
  requireNamespace("fastICA")
  space <- fastICA::fastICA(t(scale(t(x))), ndim)$S
  process_dimred(space, rownames(x))
}

dimred_lle <- function(x, ndim = 3) {
  requireNamespace("lle")
  poss_k <- lle::calc_k(t(scale(t(x))), ndim)
  k <- poss_k$k[which.min(poss_k$rho)]
  space <- lle::lle(t(scale(t(x))), ndim, k)$Y
  process_dimred(space, rownames(x))
}
