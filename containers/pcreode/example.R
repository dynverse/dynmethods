set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/pcreode",
  num_cells = 199,
  num_features = 101,
  model = "linear"
)
params <- list()
