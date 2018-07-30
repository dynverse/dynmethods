library(tidyverse)
library(furrr)

plan(multiprocess)

definition_files <- list.files("containers", pattern = "definition.yml", recursive = TRUE, full.names = TRUE)

#' @examples
#' method <- "angle"
#' file <- paste0("containers/", method, "/definition.yml")

# rebuild all dockers in the 'containers' folder
future_map(definition_files, function(file) {
  definition <- yaml::read_yaml(file)
  docker_repo <- definition$docker_repository
  folder <- stringr::str_replace(file, "/definition.yml$", "")

  # processx::run("docker", args = c("pull", docker_repo), echo = TRUE)
  processx::run("docker", args = c("build", folder, "-t", docker_repo), echo = TRUE)
  processx::run("docker", args = c("push", docker_repo), echo = TRUE)
})

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_features = 300, model = "binary_tree")
#' traj <- dynwrap::infer_trajectory(data, method, verbose = TRUE)
#' dynplot::plot_graph(traj)
