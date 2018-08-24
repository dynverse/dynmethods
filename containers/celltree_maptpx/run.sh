library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

requireNamespace("cellTree")

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "binary_tree") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/celltree_gibbs/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)
#' params$tot_iter <- 30 # override for testing

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

expression <- data$expression

start_cell <-
  if (!is.null(data$start_id)) {
    sample(data$start_id, 1)
  } else {
    NULL
  }

if (params$rooting_method == "null") {
  params$rooting_method <- NULL
}

if (is.null(params$num_topics)) {
  num_topics <- seq(params$num_topics_lower, params$num_topics_upper)
} else {
  num_topics <- params$num_topics
}

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# infer the LDA model
lda_out <- cellTree::compute.lda(
  t(expression) + min(expression) + 1,
  k.topics = num_topics,
  method = params$method,
  log.scale = FALSE,
  sd.filter = params$sd_filter,
  tot.iter = params$tot_iter,
  tol = params$tolerance
)

# check whether there is prior information available
start.group.label <- NULL
grouping <- NULL

if(!is.null(data$groups_id)) {
  grouping <-
    data$groups_id %>%
    dplyr::slice(match(cell_id, rownames(expression))) %>%
    pull(group_id)
  if(!is.null(data$start_cell)) {
    start.group.label <-
      data$groups_id %>% filter(cell_id == data$start_cell) %>%
      pull(group_id)
  }
}

# construct the backbone tree
mst_tree <- cellTree::compute.backbone.tree(
  lda.results = lda_out,
  grouping = grouping,
  start.group.label = start.group.label,
  absolute.width = params$absolute_width,
  width.scale.factor = params$width_scale_factor,
  outlier.tolerance.factor = params$outlier_tolerance_factor,
  rooting.method = params$rooting_method,
  only.mst = FALSE,
  merge.sequential.backbone = FALSE
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# simplify sample graph to just its backbone
cell_graph <- igraph::as_data_frame(mst_tree, "edges") %>%
  dplyr::select(from, to, length = weight) %>%
  mutate(
    from = rownames(expression)[from],
    to = rownames(expression)[to],
    directed = FALSE
  )
to_keep <- igraph::V(mst_tree)$is.backbone %>%
  setNames(rownames(expression))

# wrap output
output <- lst(
  cell_ids = rownames(expression),
  cell_graph,
  to_keep,
  is_directed = FALSE,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/ti/output/output.rds")
