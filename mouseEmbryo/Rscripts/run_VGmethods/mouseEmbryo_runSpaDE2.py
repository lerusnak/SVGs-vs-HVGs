"""Run SpatialDE2 on Mouse Embryo"""
import warnings
warnings.filterwarnings("ignore")

import os
absolute_dir = '/projectnb/weber-lr/'
os.chdir(absolute_dir)

import pandas as pd
import numpy as np
import argparse
import scanpy as sc
import tracemalloc
import SpatialDE as sd
from tqdm import tqdm 
from functools import partialmethod

import NaiveDE

tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

    
mouseEmbh5ad = '/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/mouseEmbryo1474343f6beb38.h5ad'

adata = sc.read_h5ad(mouseEmbh5ad)
adata

adata.obsm["spatial"] = adata.obs[['spatial0', 'spatial1']]
adata
    
# normalization - previously done in preprocessing script
#sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=[10])

#counts = sc.get.obs_df(adata, 
                      # keys=list(adata.var_names), 
                      # use_raw=False, layer='logcounts')

# total_counts = sc.get.obs_df(adata, keys=["total_counts"])
#norm_expr = NaiveDE.stabilize(counts.T).T
#adata.X = NaiveDE.regress_out(total_counts, norm_expr.T, "np.log(total_counts)").T
    
# run SpatialDE2
df_res = sd.fit(adata, normalized=True, control=None)
df_res.set_index("gene", inplace=True)
df_res = df_res.loc[adata.var_names]

# For downsream analysis 
df_res.to_csv(path_or_buf='/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/mouseEmbryo_SpaDE2_results.csv')

