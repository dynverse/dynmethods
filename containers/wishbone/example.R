set.seed(1)
data <- dyntoy::generate_dataset(
  id = "wishbone_example",
  num_cells = 99,
  num_features = 101,
  model = "bifurcating"
)
params <- list()
