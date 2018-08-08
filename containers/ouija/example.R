set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/ouija",
  num_cells = 31,
  num_features = 29,
  model = "linear"
)
params <- list(
  iter = 10
)
