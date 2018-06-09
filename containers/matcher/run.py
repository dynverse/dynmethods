import pandas as pd
import numpy as np
import json
import os

from pymatcher import matcher

# load data
expression = pd.read_csv("/input/expression.csv", index_col=[0])
params = json.load(open("/input/params.json", "r"))

# run matcher
m = matcher.MATCHER([expression.values])
m.infer(quantiles = params["quantiles"], method = params["method"])

# extract pseudotime
pseudotime = pd.DataFrame({
  "cell_id":expression.index,
  "pseudotime":m.master_time[0][:, 0]
})
pseudotime.to_csv("/output/pseudotime.csv", index=False)
