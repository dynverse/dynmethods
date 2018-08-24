import pandas as pd
import numpy as np
import json
import sys
from GPfates import GPfates

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
p = json.load(open("/ti/input/params.json", "r"))
end_n = json.load(open("/ti/input/end_n.json"))[0]
expression = pd.read_csv("/ti/input/expression.csv", index_col=0, header=0)
expression = expression[(expression > p["log_expression_cutoff"]).sum(1) >= p["min_cells_expression_cutoff"]]

if expression.shape[0] == 0:
  raise ValueError("Expression preprocessing filtered out all cells")

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
cellinfo = pd.DataFrame({"cell_id":expression.index.tolist()}, index = expression.index.tolist())
m = GPfates.GPfates(cellinfo, expression.T)

# dimensionality reduction
m.dimensionality_reduction()
m.store_dr(dims=range(p["ndim"])) # store the dr in the sample table (m.s), so it can be used in the gplvm

# infer pseudotime
m.infer_pseudotime(s_columns=["bgplvm_" + str(i) for i in range(p["ndim"])]) # use the first two components to infer pseudotime

# model different fates
m.model_fates(C=end_n)

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process and save output                                                 ####
# cell ids
cell_ids = pd.DataFrame({
  "cell_ids": m.s.pseudotime.index
})
cell_ids.to_csv("/ti/output/cell_ids.csv", index=False)

# pseudotime
pseudotime = m.s.pseudotime.reset_index()
pseudotime.columns = ["cell_id", "pseudotime"]
pseudotime.to_csv("/ti/output/pseudotime.csv", index=False)

# dimred
dimred = pd.DataFrame(m.dr_models["bgplvm"].X.mean[:,:].tolist())
dimred["cell_id"] = m.s.pseudotime.index.tolist()
dimred.to_csv("/ti/output/dimred.csv", index=False)

# progressions
end_state_probabilities = pd.DataFrame(m.fate_model.phi)
end_state_probabilities["cell_id"] = m.s.pseudotime.index.tolist()
end_state_probabilities.to_csv("/ti/output/end_state_probabilities.csv", index=False)

# timings
json.dump(checkpoints, open("/ti/output/timings.json", "w"))
