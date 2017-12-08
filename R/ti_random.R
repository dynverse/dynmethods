#' Description for random
#'
#' @importFrom dynplot plot_default
#'
#' @export
description_random <- function() create_description(
  name = "random",
  short_name = "random",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeNumericParam(id = "dummy_param", lower = 0, default = 0.5, upper = 1)
  ),
  properties = c(),
  run_fun = run_random,
  plot_fun = dynplot::plot_default
)

run_random <- function(counts, dummy_param = .5) {
  # generate network
  milestone_ids <- paste0("milestone_", seq_len(10))

  gr <- igraph::ba.game(10)
  milestone_network <- igraph::as_data_frame(gr) %>%
    mutate(
      from = paste0("milestone_", from),
      to = paste0("milestone_", to),
      length = 1,
      directed = FALSE
    )

  # put cells on random edges of network
  cell_ids <- rownames(counts)

  progressions <- data.frame(
    cell_id = cell_ids,
    milestone_network[sample.int(nrow(milestone_network), length(cell_ids), replace = TRUE), 1:2],
    percentage = runif(length(cell_ids))
  )

  # return output
  wrap_prediction_model(
    cell_ids = cell_ids,
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions
  )
}

