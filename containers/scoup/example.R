set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/scoup",
  num_cells = 49,
  num_features = 51,
  model = "tree"
)
params <- list()
