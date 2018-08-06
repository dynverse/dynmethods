set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/fateid",
  num_cells = 200,
  num_features = 101,
  model = "multifurcating"
)
params <- list(force = TRUE)
