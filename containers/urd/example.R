set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/urd",
  num_cells = 300,
  num_features = 40,
  model = "bifurcating"
)
params <- list(
  max_iter = 3
)
