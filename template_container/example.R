#!/usr/local/bin/Rscript

# generate dataset with certain seed
set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/template",
  num_cells = 200,
  num_features = 300,
  model = "tree"
)

# add method specific args (if needed)
data$params <- list()

# write example dataset to file
file <- commandArgs(trailingOnly = TRUE)[[1]]
dynutils::write_h5(data, file)
