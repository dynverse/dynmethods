library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(igraph)

library(CellTrails)

checkpoints <- list()

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 500, num_genes = 300, model = "binary_tree") %>% c(., .$prior_information)
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
tfeat1 <- filterTrajFeaturesByDL(sce, threshold = params$threshold_dl, show_plot = FALSE)
tfeat2 <- filterTrajFeaturesByCOV(sce, threshold = params$threshold_cov, show_plot = FALSE)
tfeat3 <- filterTrajFeaturesByFF(sce, threshold = params$threshold_ff, show_plot = FALSE)

trajFeatureNames(sce) <- Reduce(intersect, list(tfeat1, tfeat2, tfeat3))

# dimensionality reduction
se <- embedSamples(sce)
d <- findSpectrum(se$eigenvalues, frac = params$frac)
latentSpace(sce) <- se$components[, d]

# find states
states(sce) <- sce %>% findStates(
  min_size = params$min_size,
  min_feat = params$min_feat,
  max_pval = params$max_pval,
  min_fc = params$min_fc
)

# construct tree
sce <- connectStates(sce, l = params$l)

# fit trajectory
# this object can contain multiple trajectories (= "components"), so we have to extract information for every one of them and combine afterwards
components <- CellTrails::trajComponents(sce)


trajectories <- map(
  seq_along(components),
  function(ix) {
    if (length(components[[ix]]) > 1) {
      traj <- selectTrajectory(sce, ix)
      fitTrajectory(traj)
    } else {
      components[[ix]]
    }
  }
)


checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Process cell graph                                                      ####

cell_ids <- sampleNames(sce)
grouping <- sce$CellTrails.state %>% as.character() %>% set_names(cell_ids)
dimred <- sce@reducedDims$CellTrails

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
  grouping,
  dimred,
  cell_graph,
  to_keep,
  timings = checkpoints
)

write_rds(output, "/output/output.rds")


traj <- wrap_data(cell_ids = cell_ids) %>% add_cell_graph(cell_graph, to_keep)
