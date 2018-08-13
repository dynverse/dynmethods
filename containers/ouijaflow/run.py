import json
import pandas as pd
from ouijaflow import ouija

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
expression = pd.read_csv("/ti/input/expression.csv", index_col=[0])
p = json.load(open("/ti/input/params.json", "r"))
features_id = json.load(open("/ti/input/features_id.json"))

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
cell_ids = pd.DataFrame({
  "cell_ids": expression.index
})
cell_ids.to_csv("/ti/output/cell_ids.csv", index=False)

pseudotime = pd.DataFrame({
  "pseudotime": z,
  "cell_id": expression.index
})
pseudotime.to_csv("/ti/output/pseudotime.csv", index=False)

# timings
json.dump(checkpoints, open("/ti/output/timings.json", "w"))
