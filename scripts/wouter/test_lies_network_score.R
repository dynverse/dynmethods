networks <- list(
  linear = tibble(from=1, to=2),
  bifurcating = tibble(from=c(1, 2, 2), to=c(2, 3, 4)),
  butterfly = tibble(from=c(1, 2, 2), to=c(2, 1, 1)),
  cycle = tibble(from=1, to=1),
  tree = dyngen:::generate_random_tree() %>% select(-singular)
) %>% map(~mutate(., directed=TRUE, length=1))

network_names <- c("cycle", "cycle")
network_names <- c("linear", "cycle")
network_names <- c("cycle", "linear")
network_names <- c("tree", "linear")
network_names <- c("tree", "tree")
network_names <- c("bifurcating", "bifurcating")
network_names <- c("bifurcating", "cycle")
network_names <- c("butterfly", "linear")
network_names <- c("cycle", "linear")

network_combinations <- expand.grid(names(networks), names(networks)) %>% {split(., seq_len(nrow(.)))} %>% parallel::mclapply(mc.cores=8, function(network_names) {
  network_name1 <- as.character(network_names[[1]])
  network_name2 <- as.character(network_names[[2]])

  net1 <- networks[[network_name1]]
  net2 <- networks[[network_name2]]

  tibble(score=dyneval:::calculatie_lies_network_score(net1, net2), network_name1=network_name1, network_name2=network_name2)
}) %>% bind_rows()

network_combinations %>% ggplot() + geom_raster(aes(network_name1, network_name2, fill=score)) + viridis::scale_fill_viridis(limits=c(0, 1))
