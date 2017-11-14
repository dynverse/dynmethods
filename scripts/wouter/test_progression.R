library(tidyverse)
library(magrittr)
library(readODS)

cellinfo_known = tibble(state=c(1,1,1,2,2,2,3,3,3), progression=c(0.1, 0.5, 0.9, 0.1, 0.5, 0.9, 0.2, 0.6, 0.9)) %>% mutate(cell=1:length(state))
statenet_known = tibble(from=c(1, 1), to=c(2, 3))
statenodes_known = tibble(state=c(1, 2, 3), maxprogression=c(1, 1, 1)) %>% mutate(overallprogression=c(0, cumsum(maxprogression)[-length(maxprogression)]))
cellinfo_known$overallprogression = left_join(cellinfo_known, statenodes_known, by="state") %>% mutate(overallprogression=overallprogression+progression) %>% .$overallprogression

palettes = c("Blues", "Greens", "Oranges", "Purples", "Reds") %>% map(~colorRamp(RColorBrewer::brewer.pal(9, .)))
cellinfo_known = cellinfo_known %>% rowwise() %>% mutate(color=as.character(map2(state, progression, ~rgb(palettes[[.x]](.y), maxColorValue=255)))) %>% ungroup()

library(ggplot2)
ggplot(cellinfo_known %>% mutate(state=factor(state))) + geom_point(aes(overallprogression, state, fill=color), pch=21, colour="black", size=5) + geom_vline(aes(xintercept=overallprogression), data=statenodes_known) + scale_fill_identity()





cellinfo_known = read_ods("scripts/data/bifurcating_known.ods", sheet="progression")
statenet_known = read_ods("scripts/data/bifurcating_known.ods", sheet="statenet")
statenodes_known = read_ods("scripts/data/bifurcating_known.ods", sheet="statenodes") %>% mutate(mergedprogression=c(0, cumsum(maxprogression)[-length(maxprogression)]))

palettes = c("Blues", "Greens", "Oranges", "Purples", "Reds") %>% map(~colorRamp(RColorBrewer::brewer.pal(9, .)))
cellinfo_known = cellinfo_known %>% rowwise() %>% mutate(color=as.character(map2(state, progression, ~rgb(palettes[[.x]](.y), maxColorValue=255)))) %>% ungroup()


plot_observation = function(observation, known) {
  plotdata = left_join(observation$cellinfo, known$cellinfo %>% set_colnames(paste0("known_", colnames(known$cellinfo))), by=c("cell"="known_cell" ))
  plotdata$mergedprogression = left_join(plotdata, observation$statenodes, by="state") %>% mutate(mergedprogression=mergedprogression+progression) %>% .$mergedprogression


  ggplot(plotdata %>% mutate(state=factor(state))) + geom_point(aes(mergedprogression, state, fill=known_color), pch=21, colour="black", size=5) + geom_vline(aes(xintercept=mergedprogression), data=observation$statenodes) + scale_fill_identity()
}




get_cell_distances = function(observation) {
  snet = dyngen:::get_branchconnected_statenet(observation$statenet, observation$statenodes)
  dyngen:::get_cell_distances(observation$cellinfo, snet)
}

shuffle_progression = function(observation) {
  observation$cellinfo$cell = sample(observation$cellinfo$cell)
  observation$cellinfo = observation$cellinfo %>% arrange(cell)
  observation
}
randomize_progression = function(observation) {
  observation$cellinfo$progression = runif(nrow(observation$cellinfo), 0, 1)
  observation$cellinfo$state = sample.int(nrow(observation$statenodes), nrow(observation$cellinfo), T)
  observation
}
randomize_sample_progression = function(observation, perc=0.1) {
  sample = sample.int(nrow(observation$cellinfo), nrow(observation$cellinfo)*perc)
  observation$cellinfo$progression[sample] = runif(length(sample), 0, 1)
  observation$cellinfo$state[sample] = sample.int(nrow(observation$statenodes), length(sample), T)
  observation
}

score = function(observation, known) {
  cell_distances_known = get_cell_distances(known)
  cell_distances_observed = get_cell_distances(observation)
  tibble(
    difference=mean(abs(cell_distances_known - cell_distances_observed)),
    cor=cor(get_cell_distances(known) %>% as.numeric, get_cell_distances(observation) %>% as.numeric)
  )
}
known = dambiutils::named_list(cellinfo=cellinfo_known, statenet=statenet_known, statenodes=statenodes_known)

cases = list(known=known, randomized=randomize_progression(known), shuffled=shuffle_progression(known), randomized_half=randomize_sample_progression(known, 0.5))
map2(cases, names(cases), ~score(.x, known) %>% mutate(case=.y)) %>% bind_rows()
map2(cases, names(cases), ~plot_observation(.x, known) + ggtitle(.y)) %>% cowplot::plot_grid(plotlist=.)



cases2 = tibble(perc=rep(seq(0, 1, 0.05), 10)) %>% rowwise() %>% mutate(case=map(perc, ~randomize_sample_progression(known, .)))
scores = map(cases2$case, ~score(., known)) %>% bind_rows()
cases2 = bind_cols(cases2, scores)
ggplot(cases2) + geom_boxplot(aes(factor(perc), cor, group=perc))

