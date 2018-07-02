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

# end_state_probabilities
end_state_probabilities <- fb$probs %>% as.data.frame() %>% rownames_to_column("cell_id")
write_csv(end_state_probabilities, "/output/end_state_probabilities.csv")

# pseudotime
pseudotimes <- map2_dfr(names(pr$trc), pr$trc, function(curve_id, trc) {
    tibble(
      cell_id = trc,
      pseudotime = seq_along(trc)/length(trc),
      curve_id = curve_id
    )
  }) %>%
  arrange(pseudotime) %>%
  group_by(cell_id) %>%
  filter(pseudotime == max(pseudotime)) %>%
  filter(row_number() == 1)

# extract dimred
dimred <- dr[[1]][[1]] %>% as.data.frame() %>% mutate(cell_id = rownames(expression))

# save to files
write_csv(pseudotimes, "/output/pseudotime.csv")
write_csv(end_state_probabilities, "/output/end_state_probabilities.csv")
write_csv(dimred, "/output/dimred.csv")

write_json(checkpoints, "/output/timings.json")
