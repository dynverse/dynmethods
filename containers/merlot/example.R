set.seed(1)
data <- dyntoy::generate_dataset(
  id = "merlot_example",
  num_cells = 99,
  num_features = 101,
  model = "tree"
)
