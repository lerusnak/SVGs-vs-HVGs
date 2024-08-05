"""Run SpatialDE2 on Human DLPFC"""
import warnings
warnings.filterwarnings("ignore")

import os
absolute_dir = '/projectnb/weber-lr/'
os.chdir(absolute_dir)

import pandas as pd
import numpy as np
import scanpy as sc
import SpatialDE as sd
from tqdm import tqdm # to show progress bar when running code 
from functools import partialmethod

tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

    
humanDLPFCh5ad = '/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/humanDLPFC2a60203e6ef786.h5ad'

adata = sc.read_h5ad(humanDLPFCh5ad)
adata

adata.obsm["spatial"] = adata.obs[['spatial0', 'spatial1']]
adata
     
# run SpatialDE2
df_res = sd.fit(adata, normalized=True, control=None)
df_res.set_index("gene", inplace=True)
df_res = df_res.loc[adata.var_names]

# For downsream analysis 
df_res.to_csv(path_or_buf='/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/humanDLPFC_SpatialDE2_results.csv')
