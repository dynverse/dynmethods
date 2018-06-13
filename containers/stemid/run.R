library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)

library(StemID)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function (expression, clustnr = 30, bootnr = 50, metric = "pearson", 
    num_cluster_method = "sat", SE.method = "Tibs2001SEmax", 
    SE.factor = 0.25, B.gap = 50, cln = 30, FUNcluster = "kmedoids", 
    dimred_method = "tsne", outminc = 0, outlg = 2, probthr = 0.001, 
    thr_lower = -10, thr_upper = -5, outdistquant = 0.95, nmode = FALSE, 
    pdishuf = 2000, pthr = 1e-04, pethr = 1e-04) 
{
    requireNamespace("StemID")
    tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
    sc <- StemID::SCseq(data.frame(t(expression), check.names = FALSE))
    sc <- sc %>% StemID::filterdata(mintotal = 1, minexpr = 0, 
        minnumber = 0, maxexpr = Inf, downsample = TRUE, dsn = 1)
    do_gap <- num_cluster_method == "gap"
    do_sat <- num_cluster_method == "sat"
    sc <- sc %>% StemID::clustexp(clustnr = clustnr, bootnr = bootnr, 
        metric = metric, do.gap = do_gap, sat = do_sat, SE.method = SE.method, 
        SE.factor = SE.factor, B.gap = B.gap, cln = cln, FUNcluster = FUNcluster)
    sammonmap <- dimred_method == "sammon"
    initial_cmd <- dimred_method == "tsne_initcmd"
    sc <- sc %>% StemID::comptsne(sammonmap = sammonmap, initial_cmd = initial_cmd)
    sc <- sc %>% StemID::findoutliers(outminc = outminc, outlg = outlg, 
        probthr = probthr, thr = 10^(thr_lower:thr_upper), outdistquant = outdistquant)
    ltr <- StemID::Ltree(sc)
    ltr <- ltr %>% StemID::compentropy()
    ltr <- ltr %>% StemID::projcells(cthr = 0, nmode = nmode)
    ltr <- ltr %>% StemID::projback(pdishuf = pdishuf, nmode = nmode)
    ltr <- ltr %>% StemID::lineagetree(pthr = pthr, nmode = nmode)
    ltr <- ltr %>% StemID::compspantree()
    tl <- tl %>% add_timing_checkpoint("method_aftermethod")
    cluster_network <- data_frame(from = as.character(ltr@ldata$m[-1]), 
        to = as.character(ltr@trl$trl$kid), length = ltr@trl$dc[cbind(from, 
            to)], directed = FALSE)
    wrap_prediction_model(cell_ids = rownames(expression)) %>% 
        add_dimred_projection(milestone_ids = rownames(ltr@ldata$cnl %>% 
            as.matrix), milestone_network = cluster_network, 
            dimred_milestones = ltr@ldata$cnl %>% as.matrix, 
            dimred = ltr@ltcoord, milestone_assignment_cells = as.character(ltr@ldata$lp) %>% 
                setNames(rownames(expression)), num_segments_per_edge = 100, 
            col_ann = ltr@sc@fcol) %>% add_timings(timings = tl %>% 
        add_timing_checkpoint("method_afterpostproc"))
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')