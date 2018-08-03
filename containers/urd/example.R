set.seed(1)
data <- dyntoy::generate_dataset(
  id = "urd_example",
  num_cells = 300,
  num_features = 40,
  model = "bifurcating"
)
params <- list(
  max_iter = 3
)
