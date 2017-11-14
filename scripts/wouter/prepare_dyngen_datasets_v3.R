library(cowplot)
library(tidyverse)
library(dyneval)
library(monocle)
library(igraph)
library(cellTree)

select = dplyr::select
slice = dplyr::slice

datasets_info <- readRDS("../dyngen/results/datasets.rds")
datasets_info <- datasets_info %>% filter(stringr::str_detect(id, "^2017_04_25.*"))

output_root_folder <- "results/output_dyngentest/"
.datasets_location <- "../dyngen/results"


methods = list(
  celltree=list(func=trainLearner.ti.celltree, params=list(num_topics=4, width_scale_factor=1.5)),
  scorpius = list(func=trainLearner.ti.scorpius, params=list(num_dimensions=3, num_clusters=5)),
  monocle = list(func=trainLearner.ti.monocle, params=list(num_dimensions=3))
)



# dataset_num <- 10
for (dataset_num in seq_len(nrow(datasets_info))) {
  dataset_id <- datasets_info$id[[dataset_num]]
  cat("Processing ", dataset_num, "/", nrow(datasets_info), ": ", dataset_id, "\n", sep="")

  data_dir <- paste0(output_root_folder, dataset_id)
  dir.create(data_dir)
  data_file <- paste0(data_dir, "/data.RData")

  if (file.exists(data_file)) {
    cat("  Skipping; dataset already executed\n")
  } else {
    dataset <- dyngen::load_dataset(dataset_id, contents = dyngen::contents_dataset(experiment = F))

    task <- with(dataset, dyneval::wrap_ti_task_data(
      ti_type = model$modulenetname,
      name = info$id,
      expression = log2(counts+1),
      milestone_ids = gs$milestone_ids,
      milestone_network = gs$milestone_network,
      milestone_percentages = gs$milestone_percentages %>% slice(match(rownames(counts), cell_id))
    ))

    pred_outputs <- map(methods, ~do.call(.$func, c(list(.task=task, .subset=NULL), .$params)))

    ## calculate EM distances
    task_emdist <- compute_emlike_dist(task)
    pred_emdists <- map(pred_outputs, compute_emlike_dist)

    coranks <- map(pred_emdists, ~compute_coranking(task_emdist, .))

    cat("  Saving data\n")
    scores <- data.frame(
      method=names(coranks),
      bind_rows(map(coranks, "summary")),
      cor = cor(task_emdist %>% as.vector, invoke(cbind, map(pred_emdists, as.vector)))[1,]
    )

    save(task, pred_outputs, task_emdist, pred_emdists, coranks, scores,
         file = data_file)
  }
}

parallel::mclapply(seq_len(nrow(datasets_info)), function(dataset_num) {
  dataset_id <- datasets_info$id[[dataset_num]]
  cat("Processing ", dataset_num, "/", nrow(datasets_info), ": ", dataset_id, "\n", sep="")

  data_dir <- paste0(output_root_folder, dataset_id)
  dir.create(data_dir)
  data_file <- paste0(data_dir, "/data.RData")

  if (file.exists(data_file)) {
    load(data_file)

    plotdata_outputs <- map(pred_outputs, plotLearnerData.ti.default)
    plotdata_task <- plotLearnerData.ti.default(task)
    plotdata_all <- c(list(original=plotdata_task), plotdata_outputs)
    emdist_all <- c(list(original=task_emdist), pred_emdists)
    output_all <- c(list(original=task), pred_outputs)


    #plotLearner.ti.default(plotdata_task)
    #map(plotdata_outputs, plotLearner.ti.default)
    #map(plotdata_all, ~plotLearner.ti.combined(plotdata_task, .)) %>% cowplot::plot_grid(plotlist=.)


    #map(coranks, ~coRanking::imageplot(.$corank))

    plsz <- 1200
    # row 1
    for(methodid in seq_along(plotdata_all)) {
      methodname <- names(plotdata_all)[[methodid]]
      letter <- letters[[methodid]]
      plotdata <- plotdata_all[[methodid]]
      pred_output <- pred_outputs[[methodname]]
      title = paste0(methodname, " prediction")

      # prediction
      png(paste0(data_dir, "/plot1", letter, ".png"), plsz, plsz, res = 300)
      print(plotLearner.ti.default(plotdata) + labs(title = title))
      dev.off()

      # prediction from original percentages
      png(paste0(data_dir, "/plot2", letter, ".png"), plsz, plsz, res = 300)
      print(plotLearner.ti.combined(plotdata_task, plotdata) + labs(title = title))
      dev.off()

      if(methodname != "original") {
        png(paste0(data_dir, "/plot3", letter, ".png"), plsz, plsz, res = 300)
        print(get(paste0("plotLearner.ti.", methodname))(pred_output))
        dev.off()
      } else {
        png(paste0(data_dir, "/plot3", letter, ".png"), plsz, plsz, res = 300)
        plot.new()
        dev.off()
      }

      plot_emdist(output_all[[methodname]], emdist_all[[methodname]], plotdata, filename = paste0(data_dir, "/plot4", letter, ".png"), width = plsz / 300, height = plsz / 300)

      if(methodname == "original") {
        png(paste0(data_dir, "/plot5", letter, ".png"), plsz, plsz, res = 300)
        print(ggplot(scores) + geom_bar(aes(method, max_lcmc), stat = "identity"))
        dev.off()

        png(paste0(data_dir, "/plot6", letter, ".png"), plsz, plsz, res = 300)
        print(ggplot(scores) + geom_bar(aes(method, cor), stat = "identity"))
        dev.off()
      } else {
        png(paste0(data_dir, "/plot5", letter, ".png"), plsz, plsz, res = 300)
        coRanking::imageplot(coranks[[methodname]]$corank)
        dev.off()

        png(paste0(data_dir, "/plot6", letter, ".png"), plsz, plsz, res = 300)
        (ggplot() + geom_point(aes(as.vector(task_emdist), as.vector(emdist_all[[methodname]])), alpha = .1, size = .5) + labs(x = "original EM dist", y = paste0(methodname, " EM dist"))) %>% print
        dev.off()
      }
    }

    system(paste0(
      "cd ", data_dir, "\n",
      "bash << HERE\n",
      "convert plot1[abcd].png -append result_a.png\n",
      "convert plot2[abcd].png -append result_b.png\n",
      "convert plot3[abcd].png -append result_c.png\n",
      "convert plot4[abcd].png -append result_d.png\n",
      "convert plot5[abcd].png -append result_e.png\n",
      "convert plot6[abcd].png -append result_f.png\n",
      "convert result_[abcdef].png +append result.png\n",
      "HERE\n"))

    file.copy(paste0(data_dir, "/result.png"), paste0(data_dir, "_plot.png"), overwrite = T)
  }
})

dat_evals <- lapply(seq_len(nrow(datasets_info)), function(dataset_num) {
  dataset_id <- datasets_info$id[[dataset_num]]
  data_dir <- paste0(output_root_folder, dataset_id)
  dir.create(data_dir)
  data_file <- paste0(data_dir, "/data.RData")
  if (!file.exists(data_file)) {
    NULL
  } else {
    vals <- load(data_file)
    environment() %>% as.list() %>% .[vals]
  }
})

scores <- dat_evals %>% map_df(~ data.frame(dataset = .$task$name, ti_type = .$task$ti_type, .$scores))
ggplot(scores) +
  geom_path(aes(max_lcmc, cor, group = dataset, colour = ti_type)) +
  geom_point(aes(max_lcmc, cor, shape = method), size = 3)

ggplot(scores) +
  ggbeeswarm::geom_beeswarm(aes(method, cor, colour = ti_type))

#
# de <- dat_evals[[1]]
# qplot(as.vector(de$task_emdist), as.vector(de$pred_emdist1))
# qplot(as.vector(de$task_emdist), as.vector(de$pred_emdist2))
#
#
# dataset <- datasets[[10]]
#
# # expression_c <- t(SCORPIUS::quant.scale(t(log2(dataset$counts + 1))))
# # expression_r <- dataset$expression[rownames(expression_c),]
# #
# # pheatmap::pheatmap(t(SCORPIUS::quant.scale(expression_c)), cluster_cols = F)
# # pheatmap::pheatmap(t(SCORPIUS::quant.scale(expression_r)), cluster_cols = F)
#
# task <- dyngen_dataset_to_task(dataset, "freiu")
# plotdata <- plotLearnerData.ti.default(task)
# plotLearner.ti.default(plotdata)
# space_milestones <- plotdata$space_milestones
#
# expr <- task$expression
# ddr <- DDRTree::DDRTree(SCORPIUS::quant.scale(t(expr)), dimensions = 20, sigma = .001, maxIter = 20, lambda = NULL, param.gamma = 10, tol = .001)
# sample_coord <- data.frame(t(ddr$Z))
# traj_coord <- data.frame(t(ddr$Y))
# colnames(sample_coord) <- colnames(traj_coord) <- paste0("Comp", seq_len(ncol(sample_coord)))
# sample_coord$type <- task$milestone_ids[apply(task$milestone_percentages[,-1], 1, which.max)]
# tree_links <- ddr$stree %>% as.matrix %>% reshape2::melt(varnames=c("from", "to"), value.name = "value") %>% filter(value != 0)
# tree_links_coord <- tree_links %>% left_join(data.frame(from = seq_len(nrow(traj_coord)), from=traj_coord)) %>% left_join(data.frame(to = seq_len(nrow(traj_coord)), to=traj_coord))
#
# ggplot() +
#   geom_segment(aes(x = from.Comp1, xend = to.Comp1, y = from.Comp2, yend = to.Comp2), tree_links_coord) +
#   geom_point(aes(Comp1, Comp2, colour = type), sample_coord) +
#   scale_colour_manual(values = setNames(space_milestones$colour, space_milestones$id)) +
#   theme(panel.background = element_rect(fill = "#777777"))
# #
