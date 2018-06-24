library(jsonlite)
library(readr)
library(dplyr)
library(purrr)
library(tibble)

library(STEMNET)

checkpoints <- list()

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')
expression <- data$expression
end_id <- data$end_id
groups_id <- data$groups_id

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

checkpoints$method_afterpreproc <- as.numeric(Sys.time())

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

# extract dimred
dimred <- plot(output)$PlotData %>%
  rownames_to_column("cell_id") %>%
  select(cell_id, x, y) %>%
  rename(comp_1 = x, comp_2 = y)

#   ____________________________________________________________________________
#   Save output & process output                                            ####
write_csv(end_state_probabilities, "/output/end_state_probabilities.csv")
write_csv(enframe(pseudotime, "cell_id", "pseudotime"), "/output/pseudotime.csv")
write_json(checkpoints, "/output/checkpoints.json")
write_csv(dimred, "/output/dimred.csv")
