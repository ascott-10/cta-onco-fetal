# python_common_scripts/codes/gene_expression.py

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

def make_umaps(sample_id, adata, umap_dir, gene_color_map, cell_type_colors, cols_to_plot, genes_list, cell_type_label_col="predicted_cell_type", palette="glasbey", plot_always=True):

    sc.settings.autosave = False  # Saves plots to figdir automatically
    sc.settings.autoshow = False # Prevents plots from displaying inline

    os.makedirs(umap_dir, exist_ok=True)

    for col in adata.obs.columns: 
        adata.obs[col] = adata.obs[col].astype("category")

    adata = add_gene_binary_columns(adata, genes_list)

    if f"{cell_type_label_col}_colors" not in adata.uns:
        adata, gene_color_map, cell_type_colors = set_obs_colors(adata, palette, cell_type_colors, gene_color_map)

    gene_counts = {}
    for g in genes_list:
        col = f"{g}_binary"
        if col not in adata.obs: 
            continue
        vals = adata.obs[col].astype(str)
        gene_counts[g] = ((vals == f"{g}_pos").sum(), len(vals))

    for col in adata.obs.columns:
        if not (col in cols_to_plot or col.endswith("_binary")): 
            continue

        umap_out_path = os.path.join(umap_dir, f"{col}_umap_{sample_id}.png")
        
        if plot_always == False:
            if os.path.exists(umap_out_path): 
                continue

        legend_loc = "on data" if col in ["leiden"] else "right margin"
        title = f"{sample_id}\n{col}"

        if col.endswith("_binary"):
            gene = col.replace("_binary", "")
            if gene in gene_counts:
                n_pos, n_total = gene_counts[gene]
                title = f"{sample_id}\n{gene} ({n_pos}/{n_total} cells positive)"

        n_categories = adata.obs[col].nunique()
        fig, ax = plt.subplots(figsize=(5 + min(n_categories * 0.25, 4), 5) if legend_loc == "right margin" else (5, 5))

        sc.pl.umap(adata, color=col, legend_loc=legend_loc, frameon=False, ax=ax, title=title, show=False)

        #ax.set_aspect("equal", adjustable="box")
        plt.tight_layout()
        plt.savefig(umap_out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return adata, gene_color_map, cell_type_colors
            

  