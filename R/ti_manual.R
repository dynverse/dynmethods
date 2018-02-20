#' Description for manual
#'
#' @importFrom dynplot plot_default
#'
#' @export
description_manual <- function() create_description(
  name = "manual",
  short_name = "manual",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "person_id", values = c("wouters", "robrechtc", "liesbetm", "helenat"), default = "wouters"),
    makeDiscreteParam(id = "dimred_id", values = c("pca", "mds", "tsne", "dm"), default = "pca"),
    makeDiscreteParam(id = "run_i", values = c("1"), default = "1")
  ),
  properties = c(),
  run_fun = run_manual,
  plot_fun = plot_manual
)

run_manual <- function(
  counts,
  task,
  person_id = "wouters",
  dimred_id = "pca",
  run_i = 1
) {

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  manual_folder <- file.path(Sys.getenv("DYNALYSIS_PATH"), "analysis/data/derived_data/manual_ti/")

  run_id <- dynutils::pritt("{dimred_id}_{person_id}_{run_i}")

  predictions_location <- dynutils::pritt("{manual_folder}/predictions_{run_id}.rds")
  tryCatch(
    predictions <- readr::read_rds(predictions_location),
    error = function(e) {
      stop(dynutils::pritt("No predictions found at {predictions_location}."))
    }
  )

  if(!(task$id %in% predictions$task_id)) {
    stop(dynutils::pritt("No prediction of {task$id} found prediction found in predictions"))
  }

  prediction <- predictions %>% filter(task_id == task$id) %>% pull(prediction) %>% first()

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  prediction %>% add_timings_to_wrapper(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_manual <- function(prediction) {
  cluster_graph <- prediction$graph_scaled
  cluster_graph %>%
    ggraph::ggraph() +
    geom_point(aes(x, y), data=prediction$space, size=1) +
    ggraph::geom_edge_link(edge_colour="red", edge_width=4) +
    ggraph::geom_node_point(size=10, color="red") +
    ggraph::geom_node_point(size=5, color="white") +
    ggraph::theme_graph()
}
