library(jsonlite)
library(readr)
library(dplyr)
library(purrr)
library(tibble)

library(FateID)

checkpoints <- list()

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')
expression <- data$expression
end_id <- data$end_id
start_id <- data$start_id
groups_id <- data$groups_id

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

checkpoints$method_afterpreproc <- as.numeric(Sys.time())

# determine end groups
grouping <- groups_id$group_id %>% factor() %>% as.numeric() %>% set_names(groups_id$cell_id)
grouping <- grouping[rownames(expression)] # make sure order of cells is consistent
end_groups <- grouping[end_id] %>% unique()

# determine start group
start_group <- grouping[start_id %>% sample(1)] %>% unique()

# check if there are two or more end groups
if (length(end_groups) < 2) {
  stop("FateID requires at least two end cell populations, but according to the prior information there are only ", length(end_groups), " end populations!")
}

# based on https://github.com/dgrun/FateID/blob/master/vignettes/FateID.Rmd
x <- as.data.frame(t(expression))
y <- grouping
tar <- end_groups

# reclassify
if (params$reclassify) {
  rc <- reclassify(x, y, tar, clthr = params$clthr, nbfactor = params$nbfactor, q = params$q)
  y  <- rc$part
  x  <- rc$xf
}

# fate bias
fb  <- fateBias(x, y, tar, z=NULL, minnr=params$minnr, minnrh=params$minnrh, nbfactor=params$nbfactor)

# dimensionality reduction
dr  <- compdr(x, z=NULL, m=params$m, k=params$k)

# principal curves
pr  <- prcurve(y, fb, dr, k=params$k, m=params$m, trthr=params$trthr, start=start_group)

#   ____________________________________________________________________________
#   Process & save output                                                   ####

# extract trajectory from principal curves
# we generate one start milestone, and several end milestones for every possible end state
start_milestone_id <- "start"

progressions <- map2_df(names(pr$trc), pr$trc, function(end_milestone_id, cell_order) {
  tibble(cell_id = cell_order, percentage = seq_along(cell_order)/length(cell_order), from = start_milestone_id, to = end_milestone_id)
}) %>%
  group_by(cell_id) %>%
  mutate(percentage = percentage / n()) # divide here by number of trajectories to keep sum(percentage) <= 1

milestone_network <- progressions %>%
  group_by(from, to) %>%
  summarise(
    length = n(),
    directed = TRUE
  )

divergence_regions <- bind_rows(
  tibble(milestone_id = start_milestone_id, divergence_id = "divergence", is_start = TRUE),
  tibble(milestone_id = unique(milestone_network$to), divergence_id = "divergence", is_start = FALSE)
)

# extract dimred
dimred <- dr[[1]][[1]] %>% as.data.frame() %>% mutate(cell_id = rownames(expression))

# save to files
write_csv(milestone_network, "/output/milestone_network.csv")
write_csv(progressions, "/output/progressions.csv")
write_csv(divergence_regions, "/output/divergence_regions.csv")
write_json(checkpoints, "/output/timings.json")
write_csv(dimred, "/output/dimred.csv")
