set.seed(1)
data <- dyntoy::generate_dataset(
  unique_id = "comp1_example",
  num_cells = 99,
  num_features = 101,
  model = "linear"
)
