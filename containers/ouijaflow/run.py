import json
import pandas as pd
from ouijaflow import ouija

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
expression = pd.read_csv("/input/expression.csv", index_col=[0])
p = json.load(open("/input/params.json", "r"))
marker_feature_ids = json.load(open("/input/marker_feature_ids.json"))

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Create trajectory                                                       ####
expression = expression.ix[:,marker_feature_ids]

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
