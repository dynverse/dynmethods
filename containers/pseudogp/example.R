set.seed(1)
data <- dyntoy::generate_dataset(
  unique_id = "pseudogp_example",
  num_cells = 30,
  num_features = 31,
  model = "tree"
)
params <- list(
  iter = 10
)
