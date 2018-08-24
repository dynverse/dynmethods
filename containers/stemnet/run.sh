library(jsonlite)
library(readr)
library(dplyr)
library(purrr)
library(tibble)

library(STEMNET)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/ti/input/data.rds')
params <- jsonlite::read_json('/ti/input/params.json')

#' @examples
#' data <- dyntoy::generate_dataset(id = "test", num_cells = 300, num_features = 300, model = "bifurcating") %>% c(., .$prior_information)
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
grouping <- groups_id %>% select(cell_id, group_id) %>% deframe() %>% .[rownames(expression)]
end_groups <- grouping[end_id] %>% unique()

# remove end groups with only one cell
end_groups <- intersect(end_groups, names(table(grouping) %>% keep(~. > 2)))

# check if there are two or more end groups
if (length(end_groups) < 2) {
  msg <- paste0("STEMNET requires at least two end cell populations, but according to the prior information there are only ", length(end_groups), " end populations!")

  if (!identical(params$force, TRUE)) {
    stop(msg)
  }

  warning(msg, "\nForced to invent some end populations in order to at least generate a trajectory")
  poss_groups <- unique(grouping)
  if (length(poss_groups) == 1) {
    new_end_groups <- stats::kmeans(expression[grouping == poss_groups,], centers = 2)$cluster
    grouping[grouping == poss_groups] <- c(poss_groups, max(grouping) + 1)[new_end_groups]
    end_groups <- new_end_groups
  } else {
    end_groups <- sample(poss_groups, 2)
  }
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
  cell_ids = names(pseudotime),
  end_state_probabilities,
  pseudotime,
  timings = checkpoints
)
#   ____________________________________________________________________________
#   Save output & process output                                            ####
write_rds(output, "/ti/output/output.rds")
