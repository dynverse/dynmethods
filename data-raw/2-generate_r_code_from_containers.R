library(tidyverse)
library(dynwrap)
library(furrr)

source("data-raw/2a-helper_functions.R")

# here we do some metaprogramming to generate the ti_{method} functions
# plan(multiprocess)

# use all the dynmethods containers, but surely containers from other sources could be added as well
containers <- c(
  list.files("containers", pattern = "definition.yml", recursive = TRUE, full.names = TRUE) %>%
    map_chr(~ yaml::read_yaml(.)$docker_repository)
)

# iterate over the containers and generate R scripts for each of them
future_map(containers, generate_file_from_container)

# don't forget to regenerate the documentation
devtools::document()
devtools::install(dependencies = FALSE)
