set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/embeddr",
  num_cells = 99,
  num_features = 101,
  model = "linear"
)
params <- list()
