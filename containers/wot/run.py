import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import pandas as pd
import json

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
expression = pd.read_csv("/input/expression.csv", index_col=[0])
p = json.load(open("/input/params.json", "r"))
timecourse_discrete = json.load(open("/input/timecourse_discrete.json", "r"))

expression.index.name = "id"
expression.to_csv("/input/expression.txt", sep = "\t")

time_df = pd.DataFrame({
  "id": expression.index,
  "day": timecourse_discrete
})
time_df.to_csv("/input/cell_days.txt", index = False, sep = "\t")


cell_set = pd.DataFrame({
  "cell_set": expression.index,
  "description": expression.index,
  "cell_ids": expression.index
})
cell_set.to_csv("/input/cell_set.gmt", index = False, sep = "\t", header = False)


checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


os.makedirs("/tmaps")
returned_value = os.system("wot optimal_transport --matrix /input/expression.txt --cell_days /input/cell_days.txt --out /tmaps/tmaps --local_pca 1 --format txt")


# ...and now what?


checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Save output                                                             ####

pseudotime = pd.DataFrame({
  "cell_id": expression.index,
  "pseudotime": refined_pseudotimes
})
pseudotime.to_csv("/output/pseudotime.csv", index=False)

# timings
json.dump(checkpoints, open("/output/timings.json", "w"))
