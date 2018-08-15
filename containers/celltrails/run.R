library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(igraph)

library(CellTrails)

checkpoints <- list()

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/ti/input/data.rds")
params <- jsonlite::read_json("/ti/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 500, num_features = 300, model = "binary_tree") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/celltrails/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression

checkpoints$method_afterpreproc <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
# steps from the vignette https://dcellwanger.github.io/CellTrails/

sce <- SingleCellExperiment(assays = list(logcounts = t(expression)))

# filter features
trajFeatureNames(sce) <- filterTrajFeaturesByDL(sce, threshold = params$threshold_dl, show_plot = FALSE)
trajFeatureNames(sce) <- filterTrajFeaturesByCOV(sce, threshold = params$threshold_cov, show_plot = FALSE)
trajFeatureNames(sce) <- filterTrajFeaturesByFF(sce, threshold = params$threshold_ff, min_expr = params$min_expr, show_plot = FALSE)

# dimensionality reduction
se <- CellTrails::embedSamples(sce)
d <- CellTrails::findSpectrum(se$eigenvalues, frac = params$frac)
CellTrails::latentSpace(sce) <- se$components[, d]

# find states
CellTrails::states(sce) <- sce %>% CellTrails::findStates(
  min_size = params$min_size,
  min_feat = params$min_feat,
  max_pval = params$max_pval,
  min_fc = params$min_fc
)

# construct tree
sce <- CellTrails::connectStates(sce, l = params$l)

# fit trajectory
# this object can contain multiple trajectories (= "components"), so we have to extract information for every one of them and combine afterwards
components <- CellTrails::trajComponents(sce)


trajectories <- map(
  seq_along(components),
  function(ix) {
    if (length(components[[ix]]) > 1) {
      traj <- CellTrails::selectTrajectory(sce, ix)
      CellTrails::fitTrajectory(traj)
    } else {
      components[[ix]]
    }
  }
)


checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Process cell graph                                                      ####

cell_ids <- CellTrails::sampleNames(sce)
grouping <- CellTrails::states(sce) %>% as.character() %>% set_names(cell_ids)
dimred <- SingleCellExperiment::reducedDim(sce, type = "CellTrails")

cell_graph <- map_dfr(
  trajectories,
  function(traj) {
    if (is.character(traj)) {
      cell_ids <- colnames(sce)[which(states(sce) == traj)]
      data_frame(
        from = cell_ids[-length(cell_ids)],
        to = cell_ids[-1],
        length = 0,
        directed = FALSE
      )
    } else {
      graph <- CellTrails:::.trajGraph(traj)
      cell_ids_graph <- igraph::vertex.attributes(graph)$sampleName
      cell_graph <- graph %>%
        igraph::as_data_frame() %>%
        mutate(
          from = cell_ids_graph[as.numeric(from)],
          to = cell_ids_graph[as.numeric(to)],
          directed = FALSE
        ) %>%
        dplyr::rename(
          length = weight
        )
    }
  }
)

to_keep <- unique(c(cell_graph$from, cell_graph$to))

output <- lst(
  cell_ids = to_keep,
  grouping,
  dimred,
  cell_graph,
  to_keep,
  timings = checkpoints
)

write_rds(output, "/ti/output/output.rds")
