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

expression = pd.read_csv("/ti/input/expression.csv", index_col=[0])
p = json.load(open("/ti/input/params.json", "r"))

expression.T.to_csv("/ti/input/expression.tsv", sep="\t")

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
cell_IDs, data, markers, cell_stages, data_tag, output_directory = PySCUBA.Preprocessing.RNASeq_preprocess(
  "/ti/input/expression.tsv",
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
  output_directory="/ti/workspace")

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process & save output                                                   ####
cell_ids = pd.DataFrame({
  "cell_ids": expression.index
})
cell_ids.to_csv("/ti/output/cell_ids.csv", index=False)

# grouping
grouping = pd.DataFrame({
  "cell_id": expression.index,
  "group_id": cluster_indices.astype(str)
})
grouping.to_csv("/ti/output/grouping.csv", index=False)

# milestone_network
milestone_network = pd.DataFrame([{"from":str(i), "to":str(j)} for i, js in parent_clusters.items() for j in js])
milestone_network["length"] = 1
milestone_network["directed"] = True
milestone_network.to_csv("/ti/output/milestone_network.csv", index=False)

# timings
json.dump(checkpoints, open("/ti/output/timings.json", "w"))
