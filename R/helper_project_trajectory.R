#' @importFrom pdist pdist
project_cells_to_segments <- function(
  cluster_network,
  cluster_space,
  sample_space,
  sample_cluster = NULL,
  num_segments_per_edge = 100,
  milestone_rename_fun = NULL
) {
  # collect information on cells
  space_df <- sample_space %>%
    as.data.frame %>%
    rownames_to_column("cell_id")

  if (!is.null(sample_cluster)) {
    space_df$label <- sample_cluster
  }

  # collect information on clusters
  centers_df <- cluster_space %>%
    as.data.frame %>%
    rownames_to_column("clus_id")

  # collect information on edges
  edge_df <- cluster_network %>%
    left_join(centers_df %>% rename(from = clus_id) %>% rename_if(is.numeric, ~ paste0("from.", .)), by = "from") %>%
    left_join(centers_df %>% rename(to = clus_id) %>% rename_if(is.numeric, ~ paste0("to.", .)), by = "to")

  # construct segments
  segment_df <- edge_df %>%
    rowwise() %>%
    do(data.frame(
      from = .$from,
      to = .$to,
      percentage = seq(0, 1, length.out = num_segments_per_edge),
      sapply(colnames(sample_space), function(x) {
        seq(.[[paste0("from.", x)]], .[[paste0("to.", x)]], length.out = num_segments_per_edge)
      }),
      stringsAsFactors = FALSE
    )) %>%
    ungroup()

  # calculate shortest segment piece for each cell
  segment_ix <- sapply(seq_len(nrow(sample_space)), function(i) {
    x <- sample_space[i,]

    # limit possible edges based on sample cluster
    if (!is.null(sample_cluster)) {
      la <- space_df$label[[i]]
      ix <- which(segment_df$from == la | segment_df$to == la)
    } else {
      ix <- seq_len(nrow(segment_df))
    }

    dis <- pdist::pdist(x, segment_df[ix,colnames(sample_space)])
    wm <- which.min(as.matrix(dis)[1,])
    ix[wm]
  })

  # construct progressions
  progressions <- data.frame(
    cell_id = rownames(sample_space),
    segment_df[segment_ix,] %>% select(from, to, percentage),
    stringsAsFactors = TRUE
  )

  # collect milestone network and ids
  milestone_network <- edge_df %>%
    select(from, to, length, directed)
  milestone_ids <- rownames(cluster_space)

  # rename milestones
  if (!is.null(milestone_rename_fun)) {
    progressions <- progressions %>% mutate_at(c("from", "to"), milestone_rename_fun)
    milestone_network <- milestone_network %>% mutate_at(c("from", "to"), milestone_rename_fun)
    milestone_ids <- milestone_rename_fun(milestone_ids)
    edge_df <- edge_df %>% mutate_at(c("from", "to"), milestone_rename_fun)
    centers_df <- centers_df %>% mutate(clus_id = milestone_rename_fun(clus_id))
    space_df <- space_df %>% mutate(label = milestone_rename_fun(label))
  }

  lst(
    milestone_ids,
    milestone_network,
    progressions,
    space_df,
    centers_df,
    edge_df
  )
}
