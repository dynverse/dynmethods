library(dynmethods)
library(dynutils)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
library(tibble)
library(ggplot2)

# # expect to fail
# out <- execute_evaluation(tasks = dyntoy::toy_tasks[5,], method = description_celltree_maptpx(), parameters = list(), timeout = 4, metrics = "auc_R_nx")
# attr(out, "extras")$.summary

# expect to run!
# method <- description_celltree_gibbs()
# out <- execute_evaluation(tasks = dyntoy::toy_tasks[5,], method = method, parameters = list(), timeout = 240, metrics = "auc_R_nx", output_model = T)
# attr(out, "extras")$.summary
# prediction <- attr(out, "extras")$.model[[1]]
# method$plot_fun(prediction)

dataset <- dynutils::extract_row_to_list(dyntoy::toy_tasks, 5)
counts <- dataset$counts
# start_cell <- sample(dataset$prior_information$start_cells, 1)
# end_cells <- dataset$prior_information$end_cells
# grouping_assignment <- dataset$prior_information$grouping_assignment


# celltree gibbs

# celltree maptpx
start_cells = NULL
grouping_assignment = NULL
method = "maptpx"
num_topics_lower = 2
num_topics_upper = 15
num_topics = num_topics_lower:num_topics_upper
sd_filter = .5
tot_iter = 1e6
tolerance = .05
width_scale_factor = 1.5

# dpt
start_cells = NULL
sigma = "local"
distance = "euclidean"
n_eigs = 20
density_norm = TRUE
n_local_lower = 5
n_local_upper = 7
w_width = .1

# embeddr
kernel = "nn"
metric = "correlation"
nn_pct = 1
eps = 1
t = 1
symmetrize = "mean"
measure_type = "unorm"
p = 2
thresh = .001
maxit = 10
stretch = 2
smoother = "smooth.spline"

# gpfates
nfates = 1
ndims = 2
log_expression_cutoff = 2
min_cells_expression_cutoff = 2
num_cores = 1
verbose = FALSE

# mfa
b = 2
iter=2000
thin=1
zero_inflation=FALSE
pc_initialise=1
prop_collapse=0
scale_input=TRUE

# monocle
reduction_method <- "DDRTree"
max_components <- 2
norm_method <- "vstExprs"
auto_param_selection <- TRUE
num_paths <- NULL

# monocle 1
reduction_method <- "ICA"
max_components <- 2
norm_method <- "vstExprs"
auto_param_selection <- TRUE
num_paths <- 1

# mpath
distMethod = "euclidean"
method = "kmeans"
numcluster = 11
diversity_cut = .6
size_cut = .05

# ouija
iter = 20
response_type = "switch"
inference_type = "hmc"
normalise_expression = TRUE

# phenopath
thin = 40
z_init = 1
model_mu = FALSE
scale_y = TRUE

# pseudogp
dimreds = names(dynmethods:::list_dimred_methods()) == "pca"
chains = 1
iter = 1000
smoothing_alpha = 10
smoothing_beta = 3
pseudotime_mean = 0.5
pseudotime_var = 1
initialise_from = "random"

# SCORPIUS
ndim = 3
k = 4
distance_method = "spearman"
thresh = .001
maxit = 10
stretch = 0
smoother = "smooth.spline"

# SCOUP
ndim = 3
nbranch = 4
max_ite1 = 10
max_ite2 = 40
alpha_min = .1
alpha_max = 100
t_min = .001
t_max = 2
sigma_squared_min = .1
thresh = .01
verbose = TRUE

# SCUBA
rigorous_gap_stats = TRUE
N_dim = 2
low_gene_threshold = 1
low_gene_fraction_max = 0.7
min_split = 15
min_percentage_split = 0.25

# SLICE
grouping_assignment = NULL
lm.method = "clustering"
model.type = "tree"
ss.method = "all"
ss.threshold = 0.25
community.method = "louvain"
cluster.method = "kmeans"
k = 0
k.max = 10
B = 100
k.opt.method = "firstmax"

# SLICER
end_cells = NULL
kmin = 10
m = 2
min_branch_len = 5
min_representative_percentage = 0.8
max_same_milestone_distance = 0.1

# slingshot
start_cell = NULL
end_cells = NULL
ndim = 3
nclus = 5
dimred_name = "pca"
shrink=1
reweight=TRUE
drop.multi=TRUE
thresh=0.001
maxit=15
stretch=2
smoother = "smooth.spline"
shrink.method = "cosine"

# stemid
clustnr = 30
bootnr = 50
metric = "pearson"
num_cluster_method = "sat"
SE.method = "Tibs2001SEmax"
SE.factor = .25
B.gap = 50
cln = 0
FUNcluster = "kmedoids"
dimred_method = "tsne"
outminc = 5
outlg = 2
probthr = 1e-3
thr_lower = -40
thr_upper = -1
outdistquant = .95
nmode = FALSE
pdishuf = 2000
pthr = .01
pethr = .01

# topslam
n_components = 2
n_neighbors = 10
linear_dims = 0
max_iters = 200
dimreds = rep(TRUE, 5)

# tscan
minexpr_percent = 0
minexpr_value = 0
cvcutoff = 0
clusternum_lower = 2
clusternum_upper = 9
modelNames = "VVV"

# waterfall
num_clusters = 10

# wishbone
knn = 10
n_diffusion_components = 2
n_pca_components = 15
markers = "~"
branch = TRUE
k = 15
num_waypoints = 50
normalize = TRUE
epsilon = 1
