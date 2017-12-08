#' Description for celltree maptpx
#' @export
description_celltree_maptpx <- function() abstract_celltree_description("maptpx")

#' Description for celltree gibbs
#' @export
description_celltree_gibbs <- function() abstract_celltree_description("Gibbs")

#' Description for celltree vem
#' @export
description_celltree_vem <- function() abstract_celltree_description("VEM")

abstract_celltree_description <- function(method) {
  par_set <- switch(
    method,
    maptpx = makeParamSet(
      makeDiscreteParam(id = "method", values = "maptpx", default = "maptpx"),
      makeIntegerParam(id = "num_topics_lower", lower = 2L, upper = 15L, default = 2L),
      makeIntegerParam(id = "num_topics_upper", lower = 2L, upper = 15L, default = 15L),
      makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
      makeNumericParam(id = "tot_iter", lower = log(10^4), upper = log(10^7), default = log(10^6), trafo = function(x) round(exp(x))),
      makeNumericParam(id = "tolerance", lower = log(.001), upper = log(.5), default = log(.05), trafo = exp),
      makeNumericParam(id = "width_scale_factor", lower = 1.01, default = 1.2, upper = 2),
      forbidden = quote(num_topics_lower > num_topics_upper)
    ),
    Gibbs = makeParamSet(
      makeDiscreteParam(id = "method", values = "Gibbs", default = "Gibbs"),
      makeIntegerParam(id = "num_topics", lower = 2L, default = 4L, upper = 15L),
      makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
      makeNumericParam(id = "tot_iter", lower = log(50), upper = log(500), default = log(200), trafo = function(x) round(exp(x))),
      makeNumericParam(id = "tolerance", lower = log(10^-7), upper = log(10^-3), default = log(10^-5), trafo = exp),
      makeNumericParam(id = "width_scale_factor", lower = log(.1), default = log(1.2), upper = log(100), trafo = exp)
    ),
    VEM = makeParamSet(
      makeDiscreteParam(id = "method", values = "VEM", default = "VEM"),
      makeIntegerParam(id = "num_topics", lower = 2L, default = 4L, upper = 15L),
      makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
      makeNumericParam(id = "tot_iter", lower = log(10^4), upper = log(10^7), default = log(10^6), trafo = function(x) round(exp(x))),
      makeNumericParam(id = "tolerance", lower = log(10^-7), upper = log(10^-3), default = log(10^-5), trafo = exp),
      makeNumericParam(id = "width_scale_factor", lower = log(.1), default = log(1.5), upper = log(100), trafo = exp)
    )
  )

  create_description(
    name = pritt("cellTree with {method}"),
    short_name = pritt("CT{method}"),
    package_loaded = c(),
    package_required = c("cellTree"),
    par_set = par_set,
    properties = c(),
    run_fun = run_celltree,
    plot_fun = plot_celltree
  )
}

#' @importFrom igraph degree distances get.vertex.attribute induced_subgraph
run_celltree <- function(
  # transcriptomics data
  expression,

  # prior information
  start_cells = NULL,
  grouping_assignment = NULL,

  # parameters
  method = "maptpx",
  num_topics_lower = 2,
  num_topics_upper = 15,
  num_topics = num_topics_lower:num_topics_upper,
  sd_filter = .5,
  tot_iter = 1e6,
  tolerance = .05,
  width_scale_factor = 1.5
) {
  requireNamespace("cellTree")

  start_cell <-
    if (!is.null(start_cells)) {
      sample(start_cells, 1)
    } else {
      NULL
    }

  # infer the LDA model
  lda_out <- cellTree::compute.lda(
    t(expression) + min(expression) + 1,
    k.topics = num_topics,
    method = method,
    log.scale = FALSE,
    sd.filter = sd_filter,
    tot.iter = tot_iter,
    tol = tolerance)

  # put the parameters for the backbones in a list,
  # for adding optional grouping_assignment and (if grouping is given) start group
  backbone_params <- list(
    lda.results = lda_out,
    width.scale.factor = width_scale_factor,
    only.mst = FALSE,
    merge.sequential.backbone = FALSE
  )

  # if these parameters are available, add them to the list
  if(!is.null(grouping_assignment)) {
    backbone_params$grouping <- grouping_assignment %>% slice(match(cell_id, rownames(expression))) %>% pull(group_id)
    if(!is.null(start_cell)) {
      backbone_params$start.group.label <- grouping_assignment %>% filter(cell_id == start_cell) %>% pull(group_id)
    }
  }

  # construct the backbone tree
  mst_tree <- do.call(cellTree::compute.backbone.tree, backbone_params)

  # simplify sample graph to just its backbone
  edges <- igraph::as_data_frame(mst_tree, "edges") %>%
    select(from, to, length = weight) %>%
    mutate(
      from = rownames(expression)[from],
      to = rownames(expression)[to],
      directed = FALSE
    )
  to_keep <- igraph::V(mst_tree)$is.backbone %>%
    setNames(rownames(expression))
  out <- dynutils::simplify_sample_graph(edges, to_keep, is_directed = FALSE)

  # extract data for visualisations
  tree <- cellTree:::.compute.tree.layout(mst_tree, ratio = 1)
  vertices <- igraph::as_data_frame(tree, "vertices") %>% as_data_frame()
  edges <- igraph::as_data_frame(tree, "edges") %>% as_data_frame()

  # wrap output
  wrap_prediction_model(
    trajectory_type = "tree",
    cell_ids = rownames(expression),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    vertices = vertices,
    edges = edges
  )
}

#' @importFrom ggforce geom_arc_bar
#' @importFrom grDevices rainbow
plot_celltree <- function(prediction) {
  # Based on cellTree::ct.plot.topics(prediction$mst_tree)
  vertices <- prediction$vertices
  edges <- prediction$edges

  # calculate pie sizes
  pie_df <- map_df(seq_len(nrow(vertices)), function(i) {
    pieval <- vertices$pie[[i]]
    data.frame(
      vertices[i,] %>% select(-pie),
      topic = paste0("Topic ", seq_along(pieval)),
      stringsAsFactors = FALSE
    ) %>% mutate(
      topic = factor(topic, levels = topic),
      value = pieval,
      arc = value * 2 * pi / sum(value),
      end = cumsum(arc),
      start = end - arc
    )
  })

  # obtain edge positioning
  edges_df <- data.frame(
    edges,
    from = vertices[edges$from,c("x","y")],
    to = vertices[edges$to,c("x","y")]
  )

  # get color scheme
  num_topics <- length(vertices$pie[[1]])
  ann_cols <- setNames(grDevices::rainbow(num_topics), paste0("Topic ", seq_len(num_topics)))

  # make pie graph plot
  g <- ggplot() +
    geom_segment(aes(x = from.x, xend = to.x, y = from.y, yend = to.y), edges_df) +
    ggforce::geom_arc_bar(aes(x0 = x, y0 = y, r0 = 0, r = size*2,
                              start = start, end = end, fill = topic, group = cell.name), data = pie_df, colour = NA) +
    scale_fill_manual(values = ann_cols) +
    scale_size_identity() +
    labs(fill = "Topic") +
    theme(legend.position = c(0.9, 0.125)) +
    guides(fill = guide_legend(ncol = ceiling(num_topics / 8)))
  process_dynplot(g, prediction$id)
}
