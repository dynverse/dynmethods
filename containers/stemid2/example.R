set.seed(1)
data <- dyntoy::generate_dataset(
  unique_id = "stemid2_example",
  num_cells = 199,
  num_features = 101,
  model = "tree"
)
params <- list(
  pdishuf = 50
)
