set.seed(1)
data <- dyntoy::generate_dataset(
  unique_id = "grandprix_example",
  num_cells = 99,
  num_features = 101,
  model = "multifurcating"
)
