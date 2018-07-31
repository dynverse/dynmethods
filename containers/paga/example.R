set.seed(1)
data <- dyntoy::generate_dataset(
  unique_id = "paga_example",
  num_cells = 99,
  num_features = 101,
  model = "tree"
)
