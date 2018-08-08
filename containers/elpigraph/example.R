set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/elpigraph",
  num_cells = 99,
  num_features = 101,
  model = "tree"
)
params <- list(
  NumNodes = 20L,
  MaxNumberOfIterations = 3
)
