import pandas as pd
import numpy as np
import h5py
import json

import scanpy.api as sc
import anndata

# load data
data = h5py.File("/input/data.h5", "r")
expression = pd.DataFrame(data['expression'][:].T, index=data['expression'].attrs['rownames'].astype(np.str))
data.close()

params = json.load(open("/input/params.json", "r"))

if "grouping_assignment" in data:
  grouping_assignment = data['grouping_assignment']
else:
  grouping_assignment = None

# create dataset
if grouping_assignment is not None:
  obs = grouping_assignment
  obs["louvain"] = obs.group_id.astyp("category")
  adata = anndata.AnnData(expression.values, obs)
else:
  adata = anndata.AnnData(expression.values)

# do neighbors and pca
sc.pp.neighbors(adata, n_neighbors = params["n_neighbors"])
sc.tl.pca(adata, n_comps = params["n_comps"])

# add grouping if not provided
if grouping_assignment is None:
  sc.tl.louvain(adata, resolution = params["resolution"])

# run paga
sc.tl.paga(adata)

# save output

# grouping
grouping = pd.DataFrame({"cell_id" : expression.index, "group_id":adata.obs.louvain})
grouping.reset_index(drop=True).to_feather("/output/grouping.feather")

# milestone network
milestone_network = pd.DataFrame(
  adata.uns["paga"]["connectivities_tree"].todense(),
  index=adata.obs.louvain.cat.categories,
  columns = adata.obs.louvain.cat.categories
).stack().reset_index()
milestone_network.columns = ["from", "to", "length"]
milestone_network = milestone_network.query("length > 0")
milestone_network["directed"] = False
milestone_network.reset_index(drop=True).to_feather("/output/milestone_network.feather")

# dimred
dimred = pd.DataFrame([x[0] for x in adata.obsm])
dimred.columns = ["comp_" + str(i) for i in range(dimred.shape[1])]
dimred["cell_id"] = expression.index
dimred.reset_index(drop=True).to_feather("/output/dimred.feather")

# dimred milestones
dimred["milestone_id"] = adata.obs.louvain.tolist()
dimred_milestones = dimred.groupby("milestone_id").mean().reset_index()
dimred_milestones.to_feather("/output/dimred_milestones.feather")
