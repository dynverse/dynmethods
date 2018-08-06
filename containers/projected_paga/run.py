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
data = h5py.File("/input/data.h5", "r")
counts = pd.DataFrame(data['counts'][:].T, index=data['counts'].attrs['rownames'].astype(np.str))
start_id = data['start_id'][:].astype(np.str)
data.close()

params = json.load(open("/input/params.json", "r"))

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

# grouping
grouping = pd.DataFrame({"cell_id": counts.index, "group_id": adata.obs.louvain})
grouping.reset_index(drop=True).to_feather("/output/grouping.feather")

# milestone network
milestone_network = pd.DataFrame(
  adata.uns["paga"]["connectivities_tree"].todense(),
  index=adata.obs.louvain.cat.categories,
  columns=adata.obs.louvain.cat.categories
).stack().reset_index()
milestone_network.columns = ["from", "to", "length"]
milestone_network = milestone_network.query("length > 0").reset_index(drop=True)
milestone_network["directed"] = False
milestone_network.to_feather("/output/milestone_network.feather")

# dimred
dimred = pd.DataFrame([x for x in adata.obsm['X_umap'].T]).T
dimred.columns = ["comp_" + str(i) for i in range(dimred.shape[1])]
dimred["cell_id"] = counts.index
dimred.reset_index(drop=True).to_feather("/output/dimred.feather")

# dimred milestones
dimred["milestone_id"] = adata.obs.louvain.tolist()
dimred_milestones = dimred.groupby("milestone_id").mean().reset_index()
dimred_milestones.to_feather("/output/dimred_milestones.feather")

# timings
timings = pd.Series(checkpoints)
timings.index.name = "name"
timings.name = "timings"
timings.reset_index().to_feather("/output/timings.feather")