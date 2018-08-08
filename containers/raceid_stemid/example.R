set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/raceid_stemid",
  num_cells = 199,
  num_features = 101,
  model = "tree"
)
params <- list(
  pdishuf = 50
)
