import pandas as pd
import numpy as np
import json
import os

import pandas as pd
import numpy as np
from GrandPrix import GrandPrix

import time as tm
checkpoints = {}


#   ____________________________________________________________________________
#   Load data                                                               ####
expression = pd.read_csv("/input/expression.csv", index_col=[0])
params = json.load(open("/input/params.json", "r"))

if os.path.exists("/input/time.csv"):
  time = pd.read_csv("/input/time.csv")
else:
  time = None

checkpoints["method_afterpreproc"] = tm.time()

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
# fit grandprix model, based on https://github.com/ManchesterBioinference/GrandPrix/blob/master/notebooks/Guo.ipynb
if time is not None:
  pt_np, var_np = GrandPrix.fit_model(
    data = expression.values,
    n_latent_dims = params["n_latent_dims"],
    n_inducing_points = params["n_inducing_points"],
    latent_var = params["latent_var"],
    latent_prior_mean = time.time,
    latent_prior_var = params["latent_prior_var"]
  )
else:
  pt_np, var_np = GrandPrix.fit_model(
    data = expression.values,
    n_latent_dims = params["n_latent_dims"],
    n_inducing_points = params["n_inducing_points"],
    latent_var = params["latent_var"]
  )

checkpoints["method_aftermethod"] = tm.time()

#   ____________________________________________________________________________
#   Process output                                                          ####
# process to end_state_probabilities output format
pseudotime = pd.DataFrame({
  "cell_id": expression.index,
  "pseudotime": pt_np[:, 0]
})

end_state_probabilities = pd.DataFrame({"cell_id":expression.index})
for i in range(pt_np.shape[1] - 1):
  split_id = "split_" + str(i)

  probabilities = pt_np[:, 1]
  probabilities = (probabilities - min(probabilities))/(max(probabilities) - min(probabilities))

  end_state_probabilities[split_id + "_1"] = probabilities
  end_state_probabilities[split_id + "_2"] = 1-probabilities

# save output
pseudotime.to_csv("/output/pseudotime.csv", index=False)
end_state_probabilities.to_csv("/output/end_state_probabilities.csv", index=False)

# timings
json.dump(checkpoints, open("/output/timings.json", "w"))
