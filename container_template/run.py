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

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####




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
