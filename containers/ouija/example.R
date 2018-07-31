set.seed(1)
data <- dyntoy::generate_dataset(
  unique_id = "ouija_example",
  num_cells = 31,
  num_features = 29,
  model = "linear"
)
params <- list(
  iter = 10
)
