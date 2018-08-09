set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/calista",
  num_cells = 300,
  num_features = 101,
  model = "bifurcating"
)
params <- list(
  runs = 3,
  max_iter = 5
)
