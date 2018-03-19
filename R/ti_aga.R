#' Description for aga
#' @export
description_aga <- function() abstract_aga_description("aga")

#' Description for agapt
#' @export
description_agapt <- function() abstract_aga_description("agapt")

abstract_aga_description <- function(method) {
  par_set <- makeParamSet(
    makeIntegerParam(id = "n_neighbours", lower = 1, default = 30, upper = 100),
    makeIntegerParam(id = "n_pcs", lower = 0, default = 50, upper = 100),
    makeIntegerParam(id = "n_dcs", lower = 2, default = 10, upper = 50),
    makeNumericParam(id = "resolution", lower = 0.1, default = 1, upper = 10),
    makeLogicalParam(id = "tree_based_confidence", default = TRUE)
  )

  run_fun <- switch(
    method,
    aga = run_aga,
    agapt = run_agapt
  )

  create_description(
    name = ifelse(method == "aga", "AGA", "AGA pseudotime"),
    short_name = method,
    package_loaded = c(),
    package_required = c("aga"),
    par_set = par_set,
    properties = c(),
    run_fun = run_fun,
    plot_fun = plot_aga
  )
}

## TODO: handle start cells (see below)
run_aga <- function(
  counts,
  grouping_assignment=NULL,
  start_cells=NULL,
  n_neighbours = 30,
  n_pcs = 50,
  n_dcs = 10,
  resolution = 1,
  tree_based_confidence = TRUE,
  verbose=FALSE,
  num_cores=1
) {
  requireNamespace("aga")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # RUN AGA
  # create temporary folder
  temp_folder <- tempfile()
  dir.create(temp_folder, recursive = TRUE)

  counts_file <- paste0(temp_folder, "/counts.csv")
  json_file <- paste0(temp_folder, "/params.json")

  tryCatch({
    # write counts to temporary file
    utils::write.table(t(counts), counts_file, sep=",")

    params <- as.list(environment())[methods::formalArgs(run_aga)]
    params <- params[names(params) != "counts"]

    if (!is.null(grouping_assignment)) {
      groups <- grouping_assignment %>% slice(match(rownames(counts), cell_id))
    } else {
      groups <- tibble(cell_id = rownames(counts))
    }
    groups_file <- paste0(temp_folder, "/groups.csv")
    utils::write.table(groups, groups_file, sep=",")

    # write parameters to temporary folder

    write(jsonlite::toJSON(params, auto_unbox = TRUE), json_file)

    # execute python script
    if (!is.null(num_cores)) {
      num_cores_str <- glue::glue(
        "export MKL_NUM_THREADS={num_cores};",
        "export NUMEXPR_NUM_THREADS={num_cores};",
        "export OMP_NUM_THREADS={num_cores}"
      )
    } else {
      num_cores_str <- "echo 'no cores'"
    }

    commands <- glue::glue(
      "cd {find.package('aga')}/venv",
      "source bin/activate",
      "{num_cores_str}",
      "python3 {find.package('aga')}/wrapper.py {temp_folder}",
      .sep = ";"
    )
    output <- dynutils::run_until_exit(commands)

    if (verbose) cat(output$output, "\n", sep="")

    ## load data
    obs <- read_csv(paste0(temp_folder, "/obs.csv")) %>% mutate(group_id = as.character(group_id))

    adj_ids <- c("aga_adjacency_tree_confidence", "aga_adjacency_full_confidence", "aga_adjacency_full_attachedness")
    adj <- map(
      adj_ids,
      function(adj_id) {
        adj <- scan(paste0(temp_folder, "/", adj_id, ".csv"))
        adj <- matrix(adj, nrow=sqrt(length(adj)), ncol=sqrt(length(adj)))
        rownames(adj) <- seq_len(nrow(adj))-1
        colnames(adj) <- seq_len(ncol(adj))-1
        adj %>%
          reshape2::melt(varnames=c("from", "to"), value.name=adj_id) %>%
          mutate_at(vars(from, to), as.character)
      }
    ) %>% bind_cols() %>% select(from, to, !!adj_ids)

    aga_out <- lst(obs, adj)

  }, finally = {
    # remove temporary output
    unlink(temp_folder, recursive = TRUE)
  })

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  cell_ids <- rownames(counts)

  if(is.null(start_cells)) {
    milestone_percentages <- tibble(
      cell_id = aga_out$obs$cell_id,
      milestone_id = aga_out$obs$group_id,
      percentage = 1
    )
    milestone_network <- aga_out$adj %>%
      mutate_at(vars(from, to), as.character) %>%
      filter(aga_adjacency_tree_confidence > 0) %>%
      mutate(length = 1) %>%
      mutate(directed=TRUE) %>%
      select(from, to, length, directed)

    milestone_ids <- unique(c(milestone_network$from, milestone_network$to, milestone_percentages$milestone_id))
  } else {
    # create network between branches
    branch_network <- aga_out$adj %>%
      mutate_at(vars(from, to), as.character) %>%
      filter(aga_adjacency_tree_confidence > 0) %>%
      select(from, to)

    # determine order of branches, based on location of root cell
    branch_graph <- igraph::graph_from_data_frame(branch_network)
    branch_order <- igraph::dfs(
      branch_graph,
      aga_out$obs %>% filter(cell_id == start_cells[[1]]) %>% pull(group_id)
    )$order %>%
      names()

    # now flip order of branch network if from branch is later than to branch
    branch_network <- branch_network %>%
      mutate(from_original = from, to_original = to) %>%
      mutate(flip = map2(from, to, ~diff(match(c(.x, .y), branch_order)) < 0)) %>%
      mutate(
        from = ifelse(flip, to_original, from_original),
        to = ifelse(flip, from_original, to_original),
      ) %>%
      select(from, to)

    # now create milestone network by giving each branch an edge, and adding a zero-length edge between each branch
    branch_ids <- unique(c(branch_network$from, branch_network$to, aga_out$obs$louvain_groups))

    milestone_network <- bind_rows(
      tibble(
        from=paste0(branch_ids, "_from"),
        to=paste0(branch_ids, "_to"),
        length=1,
        directed=TRUE
      ),
      branch_network %>% mutate(from = paste0(from, "_to"), to = paste0(to, "_from"), length=0, directed=TRUE)
    )

    progressions <- aga_out$obs %>%
      mutate(from = paste0(group_id, "_from"), to = paste0(group_id, "_to")) %>%
      group_by(group_id) %>%
      mutate(percentage = (aga_pseudotime - min(aga_pseudotime))/(max(aga_pseudotime) - min(aga_pseudotime))) %>%
      ungroup()

    progressions <- progressions %>%
      select(cell_id, from, to, percentage)

    milestone_ids <- unique(c(milestone_network$from, milestone_network$to, milestone_percentages$milestone_id))

    milestone_percentages <- dynwrap::convert_progressions_to_milestone_percentages(
      cell_ids,
      milestone_ids,
      milestone_network,
      progressions
    )
  }

  divergence_regions <- milestone_network %>%
    group_by(from) %>%
    filter(n() > 1) %>%
    mutate(divergence_id = from) %>%
    gather(fromto, milestone_id, from, to) %>%
    mutate(is_start = fromto == "from") %>%
    select(divergence_id, milestone_id, is_start) %>%
    distinct()

  prediction <- wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>%
    add_trajectory_to_wrapper(
      milestone_ids = milestone_ids,
      milestone_network = milestone_network,
      milestone_percentages = milestone_percentages,
      divergence_regions = divergence_regions,
      aga_out = aga_out
    )

  prediction %>%
    add_timings_to_wrapper(
      tl %>% add_timing_checkpoint("method_afterpostproc")
    )
}

run_agapt <- function(counts, start_cells) {
  invoke(run_aga, environment())
}
formals(run_agapt) <- c(
  formals(run_agapt),
  formals(run_aga)[!names(formals(run_aga)) %in% names(formals(run_agapt))]
)


plot_aga <- function(prediction) {
  requireNamespace("ggraph")
  requireNamespace("tidygraph")

  milestone_graph <- prediction$aga_out$adj %>%
    tidygraph::as_tbl_graph() %>%
    tidygraph::activate(edges) %>%
    filter(aga_adjacency_full_attachedness > 0) %>%
    # remove duplicates, retain the edge with higher weight
    arrange(-aga_adjacency_tree_confidence, -aga_adjacency_full_confidence) %>%
    mutate(edge_id = map2_chr(from, to, ~paste0(sort(c(.x, .y)), collapse="#"))) %>%
    group_by(edge_id) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    # add edge type
    mutate(
      edge_type = ifelse(aga_adjacency_tree_confidence > 0, "tree", ifelse(aga_adjacency_full_confidence > 0, "full", "attached"))
    )
  layout <- ggraph::create_layout(
    milestone_graph,
    "fr",
    weights = milestone_graph %>% tidygraph::activate(edges) %>% pull(aga_adjacency_tree_confidence) %>% {.+1}
  )
  layout %>%
    ggraph::ggraph() +
      ggraph::geom_edge_link(aes(edge_linetype = edge_type, alpha=edge_type)) +
      ggraph::geom_edge_link(aes(x=x+(xend-x)/2, y=y+(yend-y)/2, xend = x+(xend-x)/1.999, yend=y+(yend-y)/1.999), arrow=arrow(length=unit(0.1, "inches"))) +
      ggraph::scale_edge_linetype_manual(values=c(tree="solid", full="longdash", attached="dashed")) +
      ggraph::scale_edge_alpha_manual(values=c(tree=1, full=1, attached=0.5)) +
      ggraph::geom_node_point(aes(color = name), size=10) +
      geom_text(aes(x, y, label=name), size=8) +
      ggraph::theme_graph() +
      theme(legend.position = "none")
}
