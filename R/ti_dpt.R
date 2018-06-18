#' Inferring trajectories with DPT
#'
#' @inherit ti_angle description
#'
#' @inheritParams destiny::DiffusionMap
#' @inheritParams destiny::DPT
#' @param sigma Diffusion scale parameter of the Gaussian kernel. A larger sigma might be necessary if the eigenvalues can not be found because of a singularity in the matrix. Must be one of:
#'   \itemize{
#'     \item{A character vector: `"local"` (default) or `"global"`,}
#'     \item{a numeric global sigma -- a global sigma will be calculated using [destiny::find_sigmas()]}
#'     \item{or a [destiny::Sigmas-class()] object.}
#'   }
#' @param ndim Number of eigenvectors/dimensions to return
#' @param n_local_lower If sigma == 'local', the `n_local_lower`:`n_local_upper` nearest neighbor(s) determine(s) the local sigma
#' @param n_local_upper See `n_local_lower`
#' @param distance A [stats::dist()] object, or a character vector specifying which distance metric to use. Allowed measures:
#'   \itemize{
#'     \item{Euclidean distance (default),}
#'     \item{cosine distance (1-corr(c_1, c_2)), or}
#'     \item{rank correlation distance (1-corr(rank(c_1), rank(c_2)))}
#'   }
#' @importFrom reshape2 melt
#' @export
ti_dpt <- create_ti_method(
  name = "DPT",
  short_name = "dpt",
  package_loaded = c("destiny"),
  package_required = c("dynutils", "reshape2"),
  parameters = list(
    sigma = list(
      type = "discrete",
      default = "local",
      values = c("local", "global"),
      description = "Diffusion scale parameter of the Gaussian kernel. A larger sigma might be necessary if the eigenvalues can not be found because of a singularity in the matrix. Must be one of:\n\\itemize{\n\\item A character vector: \\code{\"local\"} (default) or \\code{\"global\"},\n\\item a numeric global sigma -- a global sigma will be calculated using \\code{\\link[destiny:find_sigmas]{destiny::find_sigmas()}}\n\\item or a \\code{\\link[destiny:Sigmas-class]{destiny::Sigmas-class()}} object.\n}"),

    distance = list(
      type = "discrete",
      default = "euclidean",
      values = c("euclidean", "cosine", "rankcor"),
      description = "A \\code{\\link[stats:dist]{stats::dist()}} object, or a character vector specifying which distance metric to use. Allowed measures:\n\\itemize{\n\\item Euclidean distance (default),\n\\item cosine distance (1-corr(c_1, c_2)), or\n\\item rank correlation distance (1-corr(rank(c_1), rank(c_2)))\n}"),
    ndim = list(
      type = "integer",
      default = 20L,
      upper = 100L,
      lower = 3L,
      description = "Number of eigenvectors/dimensions to return"),

    density_norm = list(
      type = "logical",
      default = TRUE,
      values = c("TRUE", "FALSE"),
      description = "logical. If TRUE, use density normalisation"),
    n_local_lower = list(
      type = "integer",
      default = 5L,
      upper = 20L,
      lower = 2L,
      description = "If sigma == 'local', the \\code{n_local_lower}:\\code{n_local_upper} nearest neighbor(s) determine(s) the local sigma"),
    n_local_upper = list(
      type = "integer",
      default = 7L,
      upper = 20L,
      lower = 2L,
      description = "See \\code{n_local_lower}"),
    w_width = list(
      type = "numeric",

      default = 0.1,
      upper = 1,
      lower = 1e-4,
      distribution = "exponential",
      rate = 1,
      description = "Window width to use for deciding the branch cutoff"),
    forbidden = "n_local_lower > n_local_upper"
  ),
  run_fun = function(
  expression,
  start_id = NULL,
  features_id = NULL,
  sigma,
  distance,
  ndim,
  density_norm,
  n_local_lower,
  n_local_upper,
  w_width
) {
  requireNamespace("destiny")

  start_cell <-
    if (!is.null(start_id)) {
      sample(start_id, 1)
    } else {
      NULL
    }

  # create n_local vector
  n_local <- seq(n_local_lower, n_local_upper, by = 1)

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run diffusion maps
  dm <- destiny::DiffusionMap(
    data = expression,
    sigma = sigma,
    distance = distance,
    n_eigs = ndim,
    density_norm = density_norm,
    n_local = n_local,
    vars = features_id
  )

  # run DPT
  dpt_params <- lst(dm, w_width)
  if (!is.null(start_cell)) {
    dpt_params$tips <- which(rownames(expression) %in% start_cell)
  }
  dpt <- do.call(destiny::DPT, dpt_params)

  # find DPT tips
  tips <- destiny::tips(dpt)
  tip_names <- rownames(expression)[tips]

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # retrieve dimred
  dimred_cells <- dpt@dm@eigenvectors %>% magrittr::set_rownames(rownames(expression)) %>% as.matrix

  # get cluster assignment
  milestone_assignment_cells <- dpt@branch[,1] %>%
    ifelse(is.na(.), 0, .) %>%
    as.character()
  branches <- sort(unique(milestone_assignment_cells))

  # calculate cluster medians
  dimred_milestones <- t(sapply(branches, function(br) colMeans(dimred_cells[milestone_assignment_cells == br,,drop=F])))

  # create star network
  milestone_network <- data_frame(
    from = "0",
    to = setdiff(branches, "0"),
    length = sqrt(rowMeans((dimred_milestones[from,] - dimred_milestones[to,])^2)),
    directed = TRUE
  )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_dimred_projection(
    milestone_ids = rownames(dimred_milestones),
    milestone_network = milestone_network,
    dimred_milestones = dimred_milestones,
    dimred = dimred_cells,
    milestone_assignment_cells = milestone_assignment_cells,
    tips = tip_names
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
},
  plot_fun = function(prediction) {
  # based on destiny::plot.DPT(prediction$dpt, col_by = "branch")

  palette <- c("#8DD3C7", "#FFED6F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#BC80BD", "#FCCDE5", "gray85", "#CCEBC5", "#FFFFB3")
  ann_cols <- c(
    setNames(c("lightgray", palette), paste0("milestone_", seq(0, length(palette)))),
    Tip = "red"
  )
  space <- prediction$dimred %>%
    data.frame() %>%
    rownames_to_column("cell_id") %>%
    mutate(
      col = ifelse(cell_id %in% prediction$tips, "Tip", prediction$milestone_assignment_cells)
    )

  g <- ggplot(space) +
    geom_point(aes(DC1, DC2, colour = col), size = 2) +
    scale_colour_manual(values = ann_cols) +
    labs(colour = "Branch") +
    theme(legend.position = c(0.9, 0.1))
  process_dynplot(g, prediction$id)
}
)
