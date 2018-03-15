#' Description for aga
#' @export
description_aga <- function() create_description(
  name = "AGA",
  short_name = "aga",
  package_loaded = c(),
  package_required = c("aga"),
  par_set = makeParamSet(
    makeNumericParam(id = "n_neighbours", lower = 1, default = 30, upper = 100),
    makeNumericParam(id = "n_pcs", lower = 0, default = 50, upper = 100),
    makeNumericParam(id = "n_dcs", lower = 2, default = 10, upper = 50)

  ),
  properties = c(),
  run_fun = run_aga,
  plot_fun = plot_aga
)

## TODO: handle start cells
run_aga <- function(
  counts,
  grouping_assignment=NULL,
  start_cells=NULL,
  n_neighbours = 30,
  n_pcs = 50,
  n_dcs = 10,
  verbose=F,
  num_cores=1
) {
  requireNamespace("aga")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  aga_args <- as.list(environment()) %>% {.[intersect(methods::formalArgs(aga::aga), names(.))]}
  aga_out <- do.call(aga::aga, aga_args)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  if(is.null(start_cells)) {
    milestone_percentages <- tibble(
      cell_id = aga_out$obs$cell_id,
      milestone_id = aga_out$obs$group_id,
      percentage = 1
    )
    milestone_network <- aga_out$adj %>%
      reshape2::melt(varnames=c("from", "to")) %>%
      mutate_at(vars(from, to), as.character) %>%
      filter(value == 1) %>%
      rename(length=value) %>%
      mutate(directed=FALSE)
    divergence_regions <- tibble(divergence_id=character(), milestone_id=character(), is_start=logical())

    milestone_ids <- unique(c(milestone_network$from, milestone_network$to, milestone_percentages$milestone_id))
    wrap_prediction_model(
      cell_ids = rownames(counts)
    ) %>%
      add_trajectory_to_wrapper(
        milestone_ids = milestone_ids,
        milestone_network = milestone_network,
        milestone_percentages = milestone_percentages,
        divergence_regions = divergence_regions
      )
  } else {
    stop("Not supported yet, have to combine pseudotimes (located in obs dataframe) with network structure. Probably will have to convert the graph to its line graph and put the cells on that by scaling the pseudotime for each branch")
  }
}


plot_aga <- function(x) {
  print("Whatever!")
}
