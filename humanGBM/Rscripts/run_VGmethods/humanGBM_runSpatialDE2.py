"""Run SpatialDE2 on Human GBM"""
import warnings
warnings.filterwarnings("ignore")

import os
absolute_dir = '/projectnb/weber-lr/'
os.chdir(absolute_dir)

import pandas as pd
import numpy as np
import scanpy as sc
import SpatialDE as sd
from tqdm import tqdm 
from functools import partialmethod

tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

    
humanGBMh5ad = '/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/humanGBM2a60207660bccf.h5ad'

adata = sc.read_h5ad(humanGBMh5ad)
adata

adata.obsm["spatial"] = adata.obs[['spatial0', 'spatial1']]
adata
      
# run SpatialDE2
df_res = sd.fit(adata, normalized=True, control=None)
df_res.set_index("gene", inplace=True)
df_res = df_res.loc[adata.var_names]

# For downsream analysis 
df_res.to_csv(path_or_buf='/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/humanGBM_SpatialDE2_results.csv')
