set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/dpt",
  num_cells = 99,
  num_features = 101,
  model = "bifurcating"
)
params <- list()
