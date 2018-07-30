library(jsonlite)
library(readr)
library(dplyr)
library(purrr)
library(tibble)

library(STEMNET)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_features = 300, model = "bifurcating") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/stemnet/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression
end_id <- data$end_id
groups_id <- data$groups_id

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# determine end groups
grouping <- groups_id %>% deframe() %>% .[rownames(expression)]
end_groups <- grouping[end_id] %>% unique()

# remove end groups with only one cell
end_groups <- intersect(end_groups, names(table(grouping) %>% keep(~.>2)))

# check if there are two or more end groups
if (length(end_groups) < 2) {
  stop("STEMNET requires at least two end cell populations, but according to the prior information there are only ", length(end_groups), " end populations with two or more cells!")
}

# create stemnet end populations
stemnet_pop <- rep(NA, nrow(expression))
stemnet_pop[which(grouping %in% end_groups)] <- grouping[which(grouping %in% end_groups)]

# run STEMNET
if (params$lambda_auto) {params$lambda <- NULL}

output <- STEMNET::runSTEMNET(
  expression,
  stemnet_pop,
  alpha = params$alpha,
  lambda = params$lambda
)

# extract pseudotime and proabilities
pseudotime <- STEMNET:::primingDegree(output)
end_state_probabilities <- output@posteriors %>% as.data.frame() %>% rownames_to_column("cell_id")

checkpoints$method_aftermethod <- as.numeric(Sys.time())

output <- lst(
  end_state_probabilities,
  pseudotime,
  timings = checkpoints
)
#   ____________________________________________________________________________
#   Save output & process output                                            ####
write_rds(output, "/output/output.rds")
