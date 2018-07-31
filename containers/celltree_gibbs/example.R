set.seed(1)
data <- dyntoy::generate_dataset(
  unique_id = "celltree_gibbs_example",
  num_cells = 99,
  num_features = 101,
  model = "tree"
)
params <- list(
  tot_iter = 10
)
