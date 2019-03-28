#!/usr/local/bin/python

import dynclipy
task = dynclipy.main()

import pandas as pd

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####

expression = task["expression"]

params = task["params"]

start_id = task["priors"]["start_id"]
if isinstance(start_id, list):
  start_id = start_id[0]

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# ... do something here ...

checkpoints["method_aftermethod"] = time.time()

# timings
timings = pd.Series(checkpoints)
timings.index.name = "name"
timings.name = "timings"
output["timings"] = timings

#   ____________________________________________________________________________
#   Save output                                                             ####

dataset = dynclipy.wrap_data(cell_ids = counts.index)
dataset.add_branch_trajectory(
  grouping = output["grouping"], 
  milestone_network = output["milestone_network"],
  branch_progressions = output["branch_progressions"],
  branches = output["branches"],
  branch_network = output["branch_network"]
)
dataset.add_dimred(
  dimred = output["dimred"]
)
dataset.write_output(task["output"])
