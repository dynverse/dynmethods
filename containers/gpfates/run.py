import pandas as pd
import numpy as np
import json
import sys
from GPfates import GPfates

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
p = json.load(open("/input/params.json", "r"))
n_end_states = json.load(open("/input/n_end_states.json"))[0]
expression = pd.read_csv("/input/expression.csv", index_col=0, header=0)
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
m.model_fates(C=n_end_states)

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process and save output                                                 ####
# pseudotime
pseudotime = m.s.pseudotime.reset_index()
pseudotime.columns = ["cell_id", "pseudotime"]
pseudotime.to_csv("/output/pseudotime.csv", index=False)

# milestone network
milestone_network = pd.DataFrame({
  "from":"M0",
  "to":["M" + str(i+1) for i in range(n_end_states)],
  "length":1,
  "directed":True
})
milestone_network.to_csv("/output/milestone_network.csv", index=False)

# dimred
dimred = pd.DataFrame(m.dr_models["bgplvm"].X.mean[:,:].tolist())
dimred["cell_id"] = m.s.pseudotime.index.tolist()
dimred.to_csv("/output/dimred.csv", index=False)

# progressions
progressions = pd.DataFrame(m.fate_model.phi)
progressions["cell_id"] = m.s.pseudotime.index.tolist()
progressions = progressions.melt(id_vars = ["cell_id"], var_name = "to", value_name = "percentage")
progressions["to"] = ["M" + str(i+1) for i in progressions["to"]]
progressions["from"] = "M0"
progressions.to_csv("/output/progressions.csv", index=False)

# divergence regions
divergence_regions = pd.DataFrame({
  "milestone_id": ["M0"] + ["M" + str(i+1) for i in range(n_end_states)],
  "divergence_id": "A",
  "is_start": [True] + [False for i in range(n_end_states)]
})
divergence_regions.to_csv("/output/divergence_regions.csv", index=False)

# timings
json.dump(checkpoints, open("/output/timings.json", "w"))
