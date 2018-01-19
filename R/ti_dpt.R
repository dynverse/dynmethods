#' Description for DPT
#' @export
description_dpt <- function() create_description(
  name = "DPT",
  short_name = "DPT",
  package_loaded = c("destiny"),
  package_required = c("dynutils", "reshape2"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "sigma", default = "local", values = c("local", "global")),
    makeDiscreteParam(id = "distance", default = "euclidean", values = c("euclidean", "cosine", "rankcor")),
    makeIntegerParam(id = "n_eigs", lower = 3L, upper = 100L, default = 20L),
    makeLogicalParam(id = "density_norm", default = TRUE),
    makeIntegerParam(id = "n_local_lower", lower = 2L, upper = 20L, default = 5L),
    makeIntegerParam(id = "n_local_upper", lower = 2L, upper = 20L, default = 7L),
    makeNumericParam(id = "w_width", lower = -4, upper = 0, default = log(.1), trafo = exp),
    forbidden = quote(n_local_lower > n_local_upper)
  ),
  properties = c(),
  run_fun = run_dpt,
  plot_fun = plot_dpt
)

#' @importFrom reshape2 melt
run_dpt <- function(expression,
                    start_cells = NULL,
                    marker_feature_ids = NULL,
                    sigma = "local",
                    distance = "euclidean",
                    n_eigs = 20,
                    density_norm = TRUE,
                    n_local_lower = 5,
                    n_local_upper = 7,
                    w_width = .1) {
  requireNamespace("destiny")

  start_cell <-
    if (!is.null(start_cells)) {
      sample(start_cells, 1)
    } else {
      NULL
    }

  # create n_local vector
  n_local <- seq(n_local_lower, n_local_upper, by = 1)

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run diffusion maps
  dm <- destiny::DiffusionMap(
    expression,
    sigma = sigma,
    distance = distance,
    n_eigs = n_eigs,
    density_norm = density_norm,
    n_local = n_local,
    vars = marker_feature_ids
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
  space <- dpt@dm@eigenvectors %>% magrittr::set_rownames(rownames(expression)) %>% as.matrix

  # get cluster assignment
  branch_assignment <- dpt@branch[,1] %>% ifelse(is.na(.), 0, .) %>% paste0("Branch", .)
  branches <- sort(unique(branch_assignment))

  # calculate cluster medians
  space_milestones <- t(sapply(branches, function(br) colMeans(space[branch_assignment == br,,drop=F])))

  # create star network
  cluster_network <- data_frame(
    from = "Branch0",
    to = setdiff(branches, "Branch0"),
    length = sqrt(rowMeans((space_milestones[from,] - space_milestones[to,])^2)),
    directed = TRUE
  )

  # TIMING: after postproc
  tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

  # transform by projecting onto cluster medians
  out <- project_cells_to_segments(
    cluster_network = cluster_network,
    cluster_space = space_milestones,
    sample_space = space,
    sample_cluster = branch_assignment
  )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    space = out$space_df,
    centers = out$centers_df,
    edge = out$edge_df,
    tips = tip_names
  ) %>% attach_timings_attribute(tl)
}

plot_dpt <- function(prediction) {
  # based on destiny::plot.DPT(prediction$dpt, col_by = "branch")

  palette <- c("#8DD3C7", "#FFED6F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#BC80BD", "#FCCDE5", "gray85", "#CCEBC5", "#FFFFB3")
  ann_cols <- c(
    setNames(c("lightgray", palette), paste0("Branch", seq(0, length(palette)))),
    Tip = "red"
  )

  g <- ggplot(prediction$space) +
    geom_point(aes(DC1, DC2, colour = ifelse(cell_id %in% prediction$tips, "Tip", label)), size = 2) +
    scale_colour_manual(values = ann_cols) +
    labs(colour = "Branch") +
    theme(legend.position = c(0.9, 0.1))
  process_dynplot(g, prediction$id)
}

