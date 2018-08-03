set.seed(1)
data <- dyntoy::generate_dataset(
  id = "cellrouter_example",
  num_cells = 200,
  num_features = 300,
  model = "tree"
)
params <- list()
