set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/wishbone",
  num_cells = 99,
  num_features = 101,
  model = "bifurcating"
)
params <- list()
