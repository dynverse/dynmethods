import pandas as pd
import json

import scimitar.models
import scimitar.morphing_mixture as mm

# load data
expression = pd.read_csv("/input/expression.csv", index_col=[0])
p = json.load(open("/input/params.json", "r"))

# run scimitar
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

# extract pseudotime
pseudotime = pd.DataFrame({
  "cell_id": expression.index,
  "pseudotime": refined_pseudotimes
})
pseudotime.to_csv("/output/pseudotime.csv", index=False)
