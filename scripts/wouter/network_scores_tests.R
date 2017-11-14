net1 <- dyngen:::generate_toy_milestone_network("linear")
net2 <- dyngen:::generate_toy_milestone_network("cycle")


### NETCOM ICS
# Does not work as expected, gives the same score for linear and cycles
net_adj1 <- net1 %>% igraph::graph_from_data_frame() %>% igraph::as_adj() %>% as.matrix()
net_adj2 <- net2 %>% igraph::graph_from_data_frame() %>% igraph::as_adj() %>% as.matrix()

alignment <- netcom::align(
  net_adj1,
  net_adj2
)

alignment$score

netcom::ics(net_adj1, net_adj2, alignment$alignment)



## NETDIST NETDIS
# Does not work for small networks (gives NaN because no large graphlets)
graphlet_size <- 2
neighbourhood_size <- 10
centred_graphlet_counts1 <- netdist::netdis_centred_graphlet_counts(net1 %>% igraph::graph_from_data_frame(), graphlet_size, neighbourhood_size)
centred_graphlet_counts2 <- netdist::netdis_centred_graphlet_counts(net2 %>% igraph::graph_from_data_frame(), graphlet_size, neighbourhood_size)

netdist::netdis(centred_graphlet_counts1, centred_graphlet_counts2, graphlet_size = graphlet_size)


net1 <- toys %>% filter(toy_category == "linear-hairy_small") %>% pull(toy) %>% first %>% .$milestone_network
net2 <- tibble(from=1, to=2)
#net2 <- dyngen:::generate_toy_milestone_network("linear")


## NETDIST NET_EMD
gdd1 <- netdist::gdd(net1 %>% igraph::graph_from_data_frame())
gdd2 <- netdist::gdd(net2 %>% igraph::graph_from_data_frame())
netdist::net_emd(gdd1, gdd2)




## OWN GA
score_map <- function(permutation, net, net_ref) {
  if(length(permutation) < nrow(net)) {
    permutation <- c(permutation, seq_len(nrow(net))[(length(permutation)+1):nrow(net)])
  }

  net_mapped <- net[permutation, permutation]

  1-sum(abs(net_mapped - net_ref))/(sum(net_mapped) + sum(net_ref))
}

complete_matrix <- function(mat, dim, fill=0) {
  mat <- rbind(mat, matrix(rep(0, ncol(mat) * (dim - nrow(mat))), ncol=ncol(mat)))
  mat <- cbind(mat, matrix(rep(0, nrow(mat) * (dim - ncol(mat))), nrow=nrow(mat)))
}

get_adjacency <- function(net, nodes=unique(c(net$from, net$to))) {
  newnet <- net %>%
    mutate(from=factor(from, levels=nodes), to=factor(to, levels=nodes)) %>%
    reshape2::acast(from~to, value.var="length", fill=0, drop=F)

  newnet[lower.tri(newnet)] = newnet[lower.tri(newnet)] + t(newnet)[lower.tri(newnet)] # make symmetric

  newnet
}


score_networks <- function(net1, net2) {
  nodes1 <- unique(c(net1$from, net1$to))
  nodes2 <- unique(c(net2$from, net2$to))

  if(length(nodes1) > length(nodes2)) {
    nodes3 <- nodes2
    net3 <- net2
    nodes2 <- nodes1
    net2 <- net1
    nodes1 <- nodes3
    net1 <- net3
  }

  net <- get_adjacency(net1)
  net_ref <- get_adjacency(net2)
  net <- complete_matrix(net, nrow(net_ref))

  # normalize weights
  net <- net/sum(net)
  net_ref <- net_ref/sum(net_ref)

  results <- GA::ga("permutation", score_map, net=net, net_ref=net_ref, min=rep(1, length(nodes1)), max=rep(length(nodes1), length(nodes1)), maxiter=10)
}

net1 <- dyngen:::generate_toy_milestone_network("linear")
net2 <- dyngen:::generate_toy_milestone_network("cycle")
net2$length <- net2$length/2

score_networks(net1, net2)
