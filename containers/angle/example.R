set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/angle",
  num_cells = 300,
  num_features = 250,
  model = "cyclic"
)
params <- list()
