import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import matplotlib
matplotlib.use('Agg')

if(not os.path.exists("graphs")):
  os.mkdir("graphs")
#%%
if(not os.path.exists("data")):
  os.mkdir("data")
#%%
if(not os.path.exists("results")):
  os.mkdir("results")

import json
import pandas as pd

print (os.getcwd())
path1=os.getcwd()

import sys
sys.path.append("/FORKS/deng_2014_python")

#%%
from forks_fcns import *

checkpoints = {}


#   ____________________________________________________________________________
#   Load data                                                               ####
counts = pd.read_csv("/ti/input/counts.csv", index_col=[0])
p = json.load(open("/ti/input/params.json", "r"))


gene_name = counts.columns
cell_names = counts.index
data = counts.values

checkpoints["method_afterpreproc"] = time()

#   ____________________________________________________________________________
#   Create trajectory                                                       ####

#%%
#pQ normalization
data1,cell_list,gene_list=normalization_pQ(data[:],p["norm_function"]=="mean",p["norm_function"]=="median",p["norm_function"]=="quantile",p["norm_quantile"])

cell_names1=cell_names.copy()
cell_names1=cell_names1[cell_list,]
#pass this way to prevent data mutation inside function: like pass by value

#%% find the number of components so as to explain the variance in the data
#pca = PCA(whiten=True)
pca = PCA()
pca.fit(data1)
mappedData=pca.fit_transform(data1)
explained_var_ratio=pca.explained_variance_ratio_
cum_sum_exp_var=np.cumsum(explained_var_ratio)
inds = np.where(cum_sum_exp_var > p["cum_sum_exp_var"])
red_dim=inds[0][0]
mappedData=mappedData[0:,0:red_dim]

#%%find number of clusters
range_clusters=range(p["min_cluster"],p["max_cluster"],1)
n_clus,max_sil,sil_scores=find_nclusters(mappedData,range_clusters)
M=range_clusters[n_clus]

mapping_type='Isomap'
init_cell=list([])

if(p["mapping_type"]=='LLE_modified'):
    n_neighbors=red_dim
else:
    n_neighbors=10

data2=data1
n_components=red_dim
isplot=False
gene_preprocessing_type='none'
mapping_params={}
mapping_params['mapping_type']=p["mapping_type"]
mapping_params['n_neighbors']=n_neighbors
mapping_params['n_components']=n_components
mapping_params['isplot']=isplot
mapping_params['gene_preprocessing_type']=gene_preprocessing_type

actual_time2=data2[:,0]

if(p["mapping_type"]=='tSNE'):
    #we first reduce the initial dimension to 50 then use tSNE on the reduced dataset
    if(data2.shape[1]>50 and red_dim<50):
        pca = PCA()
        pca.fit(data2)
        mappedData=pca.fit_transform(data2)
        mappedData=mappedData[0:,0:50]
        X_reduced=mapping(mappedData,actual_time2,mapping_params)
    elif(data2.shape[1]>2*red_dim and red_dim>50):
        pca = PCA()
        pca.fit(data2)
        mappedData=pca.fit_transform(data2)
        mappedData=mappedData[0:,0:2*red_dim]
        X_reduced=mapping(mappedData,actual_time2,mapping_params)
    else:
        X_reduced=mapping(data2,actual_time2,mapping_params)
else:
    X_reduced=mapping(data2,actual_time2,mapping_params)

print ('2 %d,%d'%X_reduced.shape)
##############################################################################
#Run steiner tree algo

if (p["initialization"]=="kmeans"):
    kmeans = KMeans(n_clusters=M, random_state=0).fit(X_reduced)
    cluster_centers_init_st=kmeans.cluster_centers_
    cluster_labels_init_st=kmeans.labels_
elif(p["initialization"]=="kmedoids"):
    kmedoids = KMedoids(n_clusters=M, random_state=0).fit(X_reduced)
    cluster_centers_init_st=kmedoids.cluster_centers_
    cluster_labels_init_st=kmedoids.labels_
else:
    np.random.seed(1)
    indices=np.random.choice(actual_time2.shape[0], M, replace=False)
    cluster_centers_init_st=X_reduced[indices,:]

print ('5 %d,%d'%X_reduced.shape)
cluster_centers_steiner,w_MST_steiner,MST_steiner,r_steiner,cluster_labels_steiner,cost_steiner,iter_steiner,all_costs_steiner=steiner_map(cluster_centers_init_st,X_reduced,p["C"],p["eta"],p["iterMax"])
print ('6 %d,%d'%data1.shape)
num_points_per_cluster_steiner=np.zeros((M,1))
for i in range(0,M):
    num_points_per_cluster_steiner[i]=np.sum(cluster_labels_steiner==i)

#find the MST connecting the steiner cluster centers
dist_mat_steiner,MST_steiner,MST_orig_steiner,w_MST_steiner=calculateMST(cluster_centers_steiner)
#search for the path ends, this step is mainly used when we do not know the starting point
end_vals_steiner,end_idxs_steiner,MST_dict_steiner=searchpathends(cluster_centers_steiner,MST_steiner)

#find the index in cluster_centers closest to the starting cell
starting_cell=np.zeros((0,))
#if you have starting cell, put labels and data as empty matrices else put starting cell as empty matrix
indices_steiner=find_closest_cluster_center(cluster_centers_steiner,X_reduced,actual_time2,starting_cell)
#search all the paths going from the cluster center nearest to starting cell
all_paths_steiner,longest_path_steiner=searchlongestpath(end_idxs_steiner,MST_dict_steiner,indices_steiner[0],num_points_per_cluster_steiner)
#perform pseudotemporal ordering from that starting cell
ordering_steiner=pseudotemporalordering(X_reduced,cluster_centers_steiner,MST_dict_steiner,all_paths_steiner,longest_path_steiner,MST_orig_steiner,cluster_labels_steiner)
#find spearman correlation with actual labels


checkpoints["method_aftermethod"] = time()

#   ____________________________________________________________________________
#   Process output & save                                                   ####
# pseudotime

cell_ids = pd.DataFrame({
  "cell_ids": cell_names1
})
cell_ids.to_csv("/ti/output/cell_ids.csv", index=False)

dimred = pd.DataFrame(
  X_reduced,
  index = cell_names1
)
dimred.index.name = "cell_id"
dimred.to_csv("/ti/output/dimred.csv", index=True)

pseudotime = pd.DataFrame({
  "pseudotime": ordering_steiner,
  "cell_id": cell_names1
})
pseudotime.to_csv("/ti/output/pseudotime.csv", index=False)

# timings
json.dump(checkpoints, open("/ti/output/timings.json", "w"))
