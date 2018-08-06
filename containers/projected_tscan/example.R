set.seed(1)
data <- dyntoy::generate_dataset(
  id = "tscan_example",
  num_cells = 99,
  num_features = 101,
  model = "tree"
)
params <- list()
