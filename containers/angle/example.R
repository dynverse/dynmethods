set.seed(1)
data <- dyntoy::generate_dataset(
  id = "angle_example",
  num_cells = 300,
  num_features = 250,
  model = "cyclic"
)
