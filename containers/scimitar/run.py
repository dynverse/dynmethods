# to avoid singularity loading the home directory python libraries
import sys
sys.path = ['/usr/local/lib/python2.7', '/usr/local/lib/python2.7/site-packages', '/usr/local/lib/python2.7/lib-old', '/usr/local/lib/python2.7/lib-dynload']

import pandas as pd
import json

import scimitar.models
import scimitar.morphing_mixture as mm

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
expression = pd.read_csv("/input/expression.csv", index_col=[0])
p = json.load(open("/input/params.json", "r"))

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
transition_model, pseudotimes = mm.morphing_gaussian_from_embedding(
  expression.values,
  fit_type='spline',
  degree=p["degree"],
  step_size=p["step_size"],
  cov_estimator=p["cov_estimator"],
  cov_reg=p["cov_reg"]
)

refined_transition_model, refined_pseudotimes = transition_model.refine(
  expression.values,
  max_iter=p["max_iter"],
  step_size=p["step_size"],
  cov_estimator=p["cov_estimator"],
  cov_reg=p["cov_reg"]
)

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
