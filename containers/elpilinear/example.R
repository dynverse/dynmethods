set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/elpilinear",
  num_cells = 99,
  num_features = 101,
  model = "linear"
)
params <- list()
