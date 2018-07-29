import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import PySCUBA
import json
import sys
import numpy as np
import os
import pandas as pd

import time
checkpoints = {}


#   ____________________________________________________________________________
#   Load data                                                               ####

expression = pd.read_csv("/input/expression.csv", index_col=[0])
p = json.load(open("/input/params.json", "r"))

# timecourse is not yet used
if os.path.exists("input/timecourse.json"):
  timecourse = json.load(open("/input/timecourse.json"))
else:
  timecourse = None

expression.T.to_csv("/input/expression.tsv", sep="\t")

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
cell_IDs, data, markers, cell_stages, data_tag, output_directory = PySCUBA.Preprocessing.RNASeq_preprocess(
  "/input/expression.tsv",
  pseudotime_mode=True,
  log_mode=False,
  N_dim=p["N_dim"],
  low_gene_threshold=p["low_gene_threshold"],
  low_gene_fraction_max=p["low_gene_fraction_max"])

centroid_coordinates, cluster_indices, parent_clusters = PySCUBA.initialize_tree(
  data,
  cell_stages,
  rigorous_gap_stats=p["rigorous_gap_stats"],
  min_split=p["min_split"],
  min_percentage_split=p["min_percentage_split"])

centroid_coordinates, cluster_indices, parent_clusters, new_tree = PySCUBA.refine_tree(
  data,
  centroid_coordinates,
  cluster_indices,
  parent_clusters,
  cell_stages,
  output_directory="/tmp")

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process & save output                                                   ####
# grouping
grouping = pd.DataFrame({
  "cell_id":expression.index,
  "group_id":cluster_indices.astype(str)
})
grouping.to_csv("/output/grouping.csv", index=False)

# milestone_network
milestone_network = pd.DataFrame([{"from":str(i), "to":str(j)} for i, js in parent_clusters.items() for j in js])
milestone_network["length"] = 1
milestone_network["directed"] = True
milestone_network.to_csv("/output/milestone_network.csv", index=False)

# extra output
pd.DataFrame(new_tree[1:], columns = new_tree[0]).to_csv("/output/new_tree.csv")

# timings
json.dump(checkpoints, open("/output/timings.json", "w"))
