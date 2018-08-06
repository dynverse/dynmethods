library(tidyverse)

rm(list = ls())

method <- "new_method"
file <- paste0("containers/", method, "/definition.yml")

# read definition and build docker
definition <- yaml::read_yaml(file)
docker_repo <- definition$docker_repository
folder <- stringr::str_replace(file, "/definition.yml$", "")
zzz <- processx::run("docker", args = c("build", folder, "-t", docker_repo), echo = TRUE)

# generate R file in dynmethods
source("data-raw/2a-helper_functions.R")
generate_file_from_container(docker_repo)
devtools::document()
devtools::install(dep = FALSE)

# rebuild if changes were made
zzz <- processx::run("docker", args = c("build", folder, "-t", docker_repo), echo = TRUE)

# try to run the method with a toy dataset
source(paste0(folder, "/example.R"))
# traj <- dynwrap::infer_trajectory(data, method, parameters = params, verbose = TRUE, debug = TRUE)
traj <- dynwrap::infer_trajectory(data, method, parameters = params, verbose = TRUE)
dynplot::plot_graph(traj)

# if it works, you can push it
# processx::run("docker", args = c("push", docker_repo), echo = TRUE)
