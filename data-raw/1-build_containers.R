library(tidyverse)
library(furrr)

plan(multiprocess)

definition_files <- list.files("containers", pattern = "definition.yml", recursive = TRUE, full.names = TRUE)

#' @examples
#' method <- "celltrails"
#' file <- paste0("containers/", method, "/definition.yml")

# rebuild all dockers in the 'containers' folder
future_map(definition_files, function(file) {
  definition <- yaml::read_yaml(file)
  docker_repo <- definition$docker_repository
  folder <- str_replace(file, "/definition.yml$", "")

  system(paste0("docker pull ", docker_repo))
  system(paste0("docker build ", folder, " -t ", docker_repo))
  system(paste0("docker push ", docker_repo))
})

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 500, num_genes = 300, model = "binary_tree")
#' traj <- dynwrap::infer_trajectory(data, method, verbose = TRUE)
#' dynplot::plot_graph(traj)
