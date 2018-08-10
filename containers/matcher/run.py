import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import pandas as pd
import numpy as np
import json
import os

from pymatcher import matcher

import time
checkpoints = {}


#   ____________________________________________________________________________
#   Load data                                                               ####
expression = pd.read_csv("/ti/input/expression.csv", index_col=[0])
params = json.load(open("/ti/input/params.json", "r"))

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
m = matcher.MATCHER([expression.values])
m.infer(quantiles = params["quantiles"], method = params["method"])

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Save output                                                             ####
cell_ids = pd.DataFrame({
  "cell_ids": expression.index
})
cell_ids.to_csv("/ti/output/cell_ids.csv", index=False)

pseudotime = pd.DataFrame({
  "cell_id":expression.index,
  "pseudotime":m.master_time[0][:, 0]
})
pseudotime.to_csv("/ti/output/pseudotime.csv", index=False)

# timings
json.dump(checkpoints, open("/ti/output/timings.json", "w"))
