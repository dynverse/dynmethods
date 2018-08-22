set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/celltrails",
  num_cells = 500,
  num_features = 300,
  model = "tree"
)
params <- list()
