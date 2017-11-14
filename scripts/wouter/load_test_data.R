library(tidyverse)
library(dyneval)

check_dependencies()

set.seed(1)
tasks <- generate_toy_datasets()
task <- dynutils::extract_row_to_list(tasks, 1)
list2env(task, .GlobalEnv)

model <- dyneval:::run_scorpius(counts)


task$geodesic_dist <- dynutils:::compute_emlike_dist(task)
model$geodesic_dist <- dynutils:::compute_emlike_dist(model)



mantel <- vegan::mantel(task$geodesic_dist, model$geodesic_dist, permutations = 1000, alternative="greater")
cor(as.numeric(task$geodesic_dist), as.numeric(model$geodesic_dist))
mantel$signif
mantel
