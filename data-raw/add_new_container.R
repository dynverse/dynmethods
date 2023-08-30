library(tidyverse)

rm(list = ls())

method <- "paga_tree"
folder <- paste0("../methods/ti_", method)

# read definition and build docker
definition <- yaml::read_yaml(paste0(folder, "/definition.yml"))
docker_repo <- definition$docker_repository
zzz <- processx::run("docker", args = c("build", folder, "-t", docker_repo), echo = TRUE)

# try to run the method with a toy dataset
source(paste0(folder, "/example.R"))

options(dynwrap_run_environment = "docker")
# traj <- dynwrap::infer_trajectory(data, docker_repo, parameters = params, verbose = TRUE, debug = TRUE)
traj <- dynwrap::infer_trajectory(data, docker_repo, parameters = params, verbose = TRUE)
dynplot::plot_graph(traj)


# # you can test whether this method can be evaluated
# eval <- dyneval::evaluate_ti_method(data, dynwrap::create_ti_methods_container(paste0("dynverse/ti_", method))(), parameters = params, metrics = c("correlation", "edge_flip", "rf_mse", "featureimp_cor"), verbose = TRUE)
# eval$summary
# dynplot::plot_graph(eval$models[[1]])
# eval$summary$error
#
# # if it works, you can push the container to docker hub
# processx::run("docker", args = c("push", docker_repo), echo = TRUE)
#
# # rebuild the singularity image
# dynbenchmark::setup_singularity_methods()
# dynwrap::pull_singularity_ti_method(docker_repo)
#
# # test the singularity image
# traj <- dynwrap::infer_trajectory(data, method, parameters = params, verbose = TRUE)
# dynplot::plot_graph(traj)
#
# # transfer it to prism
# qsub::rsync_remote(
#   remote_src = FALSE,
#   path_src = dynbenchmark::derived_file(c("dynverse/", method, ".simg"), "03-method_characterisation/singularity_images", remote = FALSE),
#   remote_dest = TRUE,
#   path_dest = dynbenchmark::derived_file(c("dynverse/", method, ".simg"), "03-method_characterisation/singularity_images", remote = TRUE),
#   verbose = TRUE
# )
