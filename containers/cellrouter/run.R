library(dplyr)
library(purrr)
library(readr)
library(feather)
library(tibble)
library(igraph)
source("/cellrouter/CellRouter_Class.R")

checkpoints <- list()

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
p <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(model = "cyclic") %>% c(., .$prior_information)
#' p <- yaml::read_yaml("containers/cellrouter/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

counts <- data$counts
grouping <- data$grouping
start_id <- data$start_id

checkpoints$method_afterpreproc <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
# steps are taken from https://github.com/edroaldo/cellrouter/blob/master/Myeloid_Progenitors/CellRouter_Paul_Tutorial.md
# but it is really hard to follow what is going on, very limited documentation
cellrouter <- CellRouter(as.data.frame(t(counts)))

cellrouter <- Normalize(cellrouter)
cellrouter <- scaleData(cellrouter)

# do pca
cellrouter <- computePCA(cellrouter, num.pcs = p$ndim_pca, seed = 42) # alarm, seed setting...!

# do tsne
if (p$max_iter == "Inf") {p$max_iter <- 100000}
cellrouter <- computeTSNE(cellrouter, num.pcs = p$ndim_tsne, seed = 42, max_iter = p$max_iter, perplexity = p$perplexity)

# louvain clustering
# cellrouter <- findClusters(cellrouter, method='model.clustering', num.pcs = 15) # this gives errors of rgba values
cellrouter <- findClusters(cellrouter, method = "graph.clustering", num.pcs = p$ndim_pca_clustering, k = p$k_clustering) # alarm, seed setting...!

# do knn
cellrouter <- buildKNN(cellrouter, k = p$k_knn, column.ann = 'population', num.pcs = p$ndim_pca_knn, sim.type = p$sim_type)

# create trajectory using start cells as source
outputdir <- "/tmp/output"
dir.create(outputdir, recursive = TRUE)
filename <- file.path(outputdir, "cell_edge_weighted_network.txt")
write.table(cellrouter@graph$edges, file = filename, sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

sources <- unique(cellrouter@sampTab$population[cellrouter@sampTab$sample_id %in% start_id])
targets <- setdiff(as.vector(cellrouter@sampTab$population), sources)

libdir <- "/cellrouter/CellRouter"
cellrouter <- findPaths(cellrouter, column='population', libdir, outputdir, method = p$distance_method_paths) # this function uses global variables...

#Preprocess trajectories
cellrouter <- processTrajectories(
  cellrouter,
  rownames(cellrouter@ndata),
  path.rank = p$ranks,
  num.cells = p$num_cells,
  neighs = p$neighs,
  column.ann = 'population',
  column.color = 'colors'
)

checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Process cell graph                                                      ####
# first get network of backbone cells
backbone_network <- map(cellrouter@paths$path, function(x) tail(head(stringr::str_split(x, "->")[[1]], -1), -1)) %>%
  map(function(order) {
    tibble(
      from = order[-length(order)],
      to = lead(order)[-length(order)],
      directed = TRUE,
      length = 1
    )
  }) %>%
  bind_rows()

# now get for every non-backbone cell the shortest backbone cell
backbone_cells <- unique(c(backbone_network$from, backbone_network$to))
nonbackbone_cells <- setdiff(rownames(counts), backbone_cells)

nonbackbone_network <- distances(cellrouter@graph$network, backbone_cells, nonbackbone_cells) %>%
  apply(2, which.min) %>%
  {backbone_cells[.]} %>%
  set_names(nonbackbone_cells) %>%
  enframe("from", "to") %>%
  mutate(length=1, directed=TRUE)

# combine to cell_graph & remove duplicated edges
cell_graph <- bind_rows(
  backbone_network,
  nonbackbone_network
) %>%
  group_by(from, to) %>%
  filter(row_number() == 1) %>%
  ungroup()

to_keep <- backbone_cells

# dimred
dimred <- cellrouter@tsne %>% as.data.frame() %>% rownames_to_column("cell_id")

# save
write_feather(cell_graph, "/output/cell_graph.feather")
write_feather(tibble(to_keep=to_keep), "/output/to_keep.feather")
write_feather(dimred, "/output/dimred.feather")

# timings
write_feather(enframe(checkpoints, "name", "timings") %>% mutate(timings = as.numeric(timings)), "/output/timings.feather")
