# python_common_scripts/codes/differential_gene_expression.py

#### Libraries ###
import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc

import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread

import matplotlib.pyplot as plt

#### Custom ####

from config import *
from codes.cell_annotation import *

def get_markers_clusters(sample_id, adata, figures_dir, tables_dir, top_n=30):

    heatmap_save_dir = os.path.join(figures_dir, "deg_analysis", "leiden", "heatmap")
    dotplot_save_dir = os.path.join(figures_dir, "deg_analysis", "leiden", "dotplot")
    markers_save_dir = os.path.join(tables_dir, "deg_analysis", "leiden")

    for d in [heatmap_save_dir, dotplot_save_dir, markers_save_dir]:
        os.makedirs(d, exist_ok=True)

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    markers_df = sc.get.rank_genes_groups_df(adata, None)
    markers_df = markers_df[~markers_df["names"].str.startswith(("RPS","RPL","MT-","MALAT1","HLA-"))]

    markers_df.to_csv(os.path.join(markers_save_dir, f"{sample_id}_deg_analysis_markers.csv"), index=False)

    top_markers = markers_df.groupby("group", observed=False).head(top_n)["names"].unique().tolist()

    top_markers = [g for g in top_markers if g in adata.var_names]

    cluster_values = adata[:, top_markers].to_df().groupby(adata.obs["leiden"]).mean().values
    cluster_values = np.nan_to_num(cluster_values)
    cluster_values = cluster_values[np.isfinite(cluster_values).all(axis=1)]

    Z = linkage(cluster_values, method="average", metric="correlation")

    sc.tl.dendrogram(adata, groupby="leiden")

    sc.pl.rank_genes_groups_dotplot(adata, show=False)
    plt.savefig(os.path.join(dotplot_save_dir, f"{sample_id}_leiden_dotplot.png"), dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.rank_genes_groups_heatmap(adata, show=False)
    plt.savefig(os.path.join(heatmap_save_dir, f"{sample_id}_leiden_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.close()

    return adata