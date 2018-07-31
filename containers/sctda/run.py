import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import matplotlib
matplotlib.use('Agg')


import pandas as pd
import json

import time
checkpoints = {}

#   ____________________________________________________________________________
#   Load data                                                               ####
expression = pd.read_csv("/input/expression.csv", index_col=[0])
p = json.load(open("/input/params.json", "r"))
timecourse_discrete = json.load(open("/input/timecourse_discrete.json", "r"))



checkpoints["method_afterpreproc"] = time.time()



#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
import scTDA

files = []
cells = []
libs = []
days = []
with open('data.txt', 'r') as f:
    for line in f:
        sp = line[:-1].split('\t')
        files.append(sp[0])
        cells.append(int(sp[1]))
        libs.append(sp[2])
        days.append(int(sp[3]))

p = scTDA.Preprocess(files, days, libs, cells, spike='ERCC-')

p.save('Embryo')
t = scTDA.TopologicalRepresentation('Embryo', lens='pca', metric='euclidean')

t.save('Embryo_dimred', 25, 0.40);

c = scTDA.UnrootedGraph('Embryo_dimred', 'Embryo.all.tsv', groups=False)

c.save(n=1000, filtercells=0, filterexp=0.0);

c = scTDA.RootedGraph('Embryo_dimred', 'Embryo.all.tsv', posgl=True, groups=False)

c.save(n=1000, filtercells=0, filterexp=0.0);

# do something with
# * Embryo_dimred.json (cluster assignment)
# * Embryo_dimred.gexf (cluster network)


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
