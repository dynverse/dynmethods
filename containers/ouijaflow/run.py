import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import json
import pandas as pd
from ouijaflow import ouija

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
expression = pd.read_csv("/input/expression.csv", index_col=[0])
p = json.load(open("/input/params.json", "r"))
features_id = json.load(open("/input/features_id.json"))

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Create trajectory                                                       ####
expression = expression.ix[:,features_id]

oui = ouija()
oui.fit(expression.values, n_iter=int(p["iter"]))

z = oui.trajectory()

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process output & save                                                   ####
# pseudotime
pseudotime = pd.DataFrame({
  "pseudotime": z,
  "cell_id": expression.index
})
pseudotime.to_csv("/output/pseudotime.csv", index=False)

# timings
json.dump(checkpoints, open("/output/timings.json", "w"))
