# avoid errors due to no $DISPLAY environment variable available when running sc.pl.paga
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import h5py
import json

import scanpy.api as sc
import anndata

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
data = h5py.File("/ti/input/data.h5", "r")
counts = pd.DataFrame(data['counts'][:].T, index=data['counts'].attrs['rownames'].astype(np.str))
start_id = data['start_id'][:].astype(np.str)
data.close()

params = json.load(open("/ti/input/params.json", "r"))

if "groups_id" in data:
  groups_id = data['groups_id']
else:
  groups_id = None

# create dataset
if groups_id is not None:
  obs = groups_id
  obs["louvain"] = obs.group_id.astyp("category")
  adata = anndata.AnnData(counts.values, obs)
else:
  adata = anndata.AnnData(counts.values)

#   ____________________________________________________________________________
#   Basic preprocessing                                                     ####

n_top_genes = min(2000, counts.shape[1])
sc.pp.recipe_zheng17(adata, n_top_genes=n_top_genes)
sc.tl.pca(adata, n_comps=params["n_comps"])
sc.pp.neighbors(adata, n_neighbors=params["n_neighbors"])

# denoise the graph by recomputing it in the first few diffusion components
if params["n_dcs"] != 0:
  sc.tl.diffmap(adata, n_comps=params["n_dcs"])

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Cluster, infer trajectory, infer pseudotime, compute dimension reduction ###

# add grouping if not provided
if groups_id is None:
  sc.tl.louvain(adata, resolution=params["resolution"])

# run paga
sc.tl.paga(adata)

# compute a layout for the paga graph
# - this simply uses a Fruchterman-Reingold layout, a tree layout or any other
#   popular graph layout is also possible
# - to obtain a clean visual representation, one can discard low-confidence edges
#   using the parameter threshold
sc.pl.paga(adata, threshold=0.01, layout='fr', show=False)

# run dpt for pseudotime information that is overlayed with paga
adata.uns['iroot'] = np.where(counts.index == start_id[0])[0][0]
sc.tl.dpt(adata)

# run umap for a dimension-reduced embedding, use the positions of the paga
# graph to initialize this embedding
if params["embedding_type"] != 'fa':
  sc.tl.draw_graph(adata, init_pos='paga')
else:
  sc.tl.umap(adata, init_pos='paga')

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process & save output                                                   ####
# cell ids
cell_ids = pd.DataFrame({
  "cell_ids": counts.index
})
cell_ids.to_feather("/ti/output/cell_ids.feather")

# grouping
grouping = pd.DataFrame({"cell_id": counts.index, "group_id": adata.obs.louvain})
grouping.reset_index(drop=True).to_feather("/ti/output/grouping.feather")

# milestone network
milestone_network = pd.DataFrame(
  adata.uns["paga"]["connectivities_tree"].todense(),
  index=adata.obs.louvain.cat.categories,
  columns=adata.obs.louvain.cat.categories
).stack().reset_index()
milestone_network.columns = ["from", "to", "length"]
milestone_network = milestone_network.query("length > 0").reset_index(drop=True)
milestone_network["directed"] = False
milestone_network.to_feather("/ti/output/milestone_network.feather")

# dimred
dimred = pd.DataFrame([x for x in adata.obsm['X_umap'].T]).T
dimred.columns = ["comp_" + str(i) for i in range(dimred.shape[1])]
dimred["cell_id"] = counts.index
dimred.reset_index(drop=True).to_feather("/ti/output/dimred.feather")

# branch progressions: the scaled dpt_pseudotime within every cluster
branch_progressions = adata.obs
branch_progressions["dpt_pseudotime"] = branch_progressions["dpt_pseudotime"].replace([np.inf, -np.inf], 1) # replace unreachable pseudotime with maximal pseudotime
branch_progressions["percentage"] = branch_progressions.groupby("louvain")["dpt_pseudotime"].apply(lambda x: (x-x.min())/(x.max() - x.min())).fillna(0.5)
branch_progressions["cell_id"] = counts.index
branch_progressions["branch_id"] = branch_progressions["louvain"].astype(np.str)
branch_progressions = branch_progressions[["cell_id", "branch_id", "percentage"]]
branch_progressions.reset_index(drop=True).to_feather("/ti/output/branch_progressions.feather")

# branches:
# - length = difference between max and min dpt_pseudotime within every cluster
# - directed = not yet correctly inferred
branches = adata.obs.groupby("louvain").apply(lambda x: x["dpt_pseudotime"].max() - x["dpt_pseudotime"].min()).reset_index()
branches.columns = ["branch_id", "length"]
branches["branch_id"] = branches["branch_id"].astype(np.str)
branches["directed"] = True
branches.to_feather("/ti/output/branches.feather")

# branch network: determine order of from and to based on difference in average pseudotime
branch_network = milestone_network[["from", "to"]]
average_pseudotime = adata.obs.groupby("louvain")["dpt_pseudotime"].mean()
for i, (branch_from, branch_to) in enumerate(zip(branch_network["from"], branch_network["to"])):
  if average_pseudotime[branch_from] > average_pseudotime[branch_to]:
    branch_network.at[i, "to"] = branch_from
    branch_network.at[i, "from"] = branch_to

branch_network.to_feather("/ti/output/branch_network.feather")

# timings
timings = pd.Series(checkpoints)
timings.index.name = "name"
timings.name = "timings"
timings.reset_index().to_feather("/ti/output/timings.feather")
