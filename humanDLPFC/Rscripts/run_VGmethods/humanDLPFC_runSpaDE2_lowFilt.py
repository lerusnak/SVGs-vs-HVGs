"""Run SpatialDE2"""
import warnings
warnings.filterwarnings("ignore")

import os
absolute_dir = '/projectnb/weber-lr/SVGs-vs-HVGs/mynewenv/lib/python3.10/site-packages'
os.chdir(absolute_dir)

import pandas as pd
import numpy as np
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import tracemalloc
import SpatialDE as sd
from tqdm import tqdm
from functools import partialmethod

import NaiveDE

tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

def main():
    args = parse_args()

    adata = sc.read_h5ad(args.input)
    
    # normalization
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=[10])

    counts = sc.get.obs_df(adata, 
                           keys=list(adata.var_names), 
                           use_raw=False, layer='counts')

    total_counts = sc.get.obs_df(adata, keys=["total_counts"])
    norm_expr = NaiveDE.stabilize(counts.T).T
    adata.X = NaiveDE.regress_out(total_counts, norm_expr.T, "np.log(total_counts)").T
    
    # run SpatialDE2
    df_res = sd.fit(adata, normalized=True, control=None)
    df_res.set_index("gene", inplace=True)
    df_res = df_res.loc[adata.var_names]

    df_res.to_csv(args.output)

if __name__ == "__main__":
    main()
