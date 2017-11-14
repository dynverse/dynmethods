library(tidyverse)
library(scent)
library(SCORPIUS)
library(org.Mm.eg.db)


expression = SCORPIUS::ginhoux$expression

#newnames = org.Mm.eg.db::org.Mm.egSYMBOL2EG %>% as.list() %>% .[colnames(expression)]
#expression = expression[, !is.na(names(newnames))]
#colnames(expression) = newnames[!is.na(names(newnames))] %>% map(~.[[1]]) %>% unlist






ppi <- read_tsv("/media/wouters/WouterSSD/downloads/BIOGRID-ORGANISM-Mus_musculus-3.4.149.mitab.txt")
ppi$from <- gsub("entrez gene/locuslink:([0-9]*)", "\\1", ppi$`#ID Interactor A`)
ppi$to <- gsub("entrez gene/locuslink:([0-9]*)", "\\1", ppi$`ID Interactor B`)
ppi$from <- mget(ppi$from, envir=org.Mm.egSYMBOL, ifnotfound=NA) %>% unlist()
ppi$to <- mget(ppi$to, envir=org.Mm.egSYMBOL, ifnotfound=NA) %>% unlist()
ppi <- ppi %>% dplyr::select(from, to) %>% filter(!is.na(from) & !is.na(to))
adj <- ppi %>% dplyr::select(from, to) %>% igraph::graph_from_data_frame(directed=FALSE) %>% igraph::as_adjacency_matrix()

scent_data <- scent::DoIntegPPI(t(expression), as.matrix(adj))

max_entropy = scent::CompMaxSR(scent_data$adjMC)
result = PRISM::qsub_lapply(colnames(scent_data$expMC), function(cell) {scent::CompSRana(scent_data$expMC[, cell], as.matrix(scent_data$adjMC), maxSR=max_entropy)}, qsub_environment = list2env(lst(max_entropy, scent_data)))


space = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation_distance(expression))
traj = SCORPIUS::infer.trajectory(space)
plot(traj$time, map_dbl(result, "sr"))

draw.trajectory.plot(space, map_dbl(result, "sr"), traj$path)
draw.trajectory.plot(space, ginhoux$sample.info$group.name)
qplot(traj$time, map_dbl(result, "sr"))
qplot(ginhoux$sample.info$group.name, map_dbl(result, "sr"))
qplot(ginhoux$sample.info$group.name, traj$time)

tibble(scorpius_time=traj$time, signaling_entropy=map_dbl(result, "sr"), group=ginhoux$sample.info$group.name) %>% GGally::ggpairs()
