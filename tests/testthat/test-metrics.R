context("Score metrics")

test_that(paste0("Check hele network score"), {
  net1 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)
  net2 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)

  expect_equal(dyneval:::calculate_lies_network_score(net1, net2), 1)

  net1 <- tibble(from=c(1, 2, 2), to=c(2, 1, 1), directed=TRUE, length=1)
  net2 <- tibble(from=c(1, 2, 2), to=c(2, 3, 4), directed=TRUE, length=1)

  expect_equal(dyneval:::calculate_lies_network_score(net1, net2), 0.5)

  net1 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)
  net2 <- tibble(from=c(1), to=c(2), directed=TRUE, length=1)

  expect_less_than(dyneval:::calculate_lies_network_score(net1, net2), 1)
})
# data("toy_tasks", package="dyntoy")
#
# for (taski in seq_len(nrow(toy_tasks))) {
#   task <- extract_row_to_list(toy_tasks, taski)
#
#   test_that(paste0("Perform dimred on trajectory with task ", task$id), {
#     g <- plot_default(task)
#     expect_is(g, "ggplot")
#
#     pdf("/dev/null")
#     print(g)
#     dev.off()
#
#     prediction <- task
#     cell_id_map <- setNames(sample(prediction$cell_ids), prediction$cell_ids)
#     prediction$milestone_percentages$cell_id <- cell_id_map[prediction$milestone_percentages$cell_id]
#     prediction$progressions$cell_id <- cell_id_map[prediction$progressions$cell_id]
#
#     g <- plot_combined(task, prediction)
#     expect_is(g, "ggplot")
#
#     pdf("/dev/null")
#     print(g)
#     dev.off()
#
#     pdf("/dev/null")
#     ph <- plot_emdist(task, task$geodesic_dist)
#     dev.off()
#     expect_is(ph$gtable, "gtable")
#   })
# }
