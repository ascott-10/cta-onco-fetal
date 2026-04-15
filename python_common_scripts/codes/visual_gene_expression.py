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
from scipy.cluster.hierarchy import linkage, leaves_list

import matplotlib.pyplot as plt
import seaborn as sns

import json

#### Custom ####

from config import *
from codes.cell_annotation import *

def make_umaps(sample_id, adata, umap_dir, gene_color_map, cell_type_colors, cols_to_plot, genes_list, cell_type_label_col="predicted_cell_type", palette="glasbey", plot_always=True):

    sc.settings.autosave = False
    sc.settings.autoshow = False

    
    os.makedirs(umap_dir, exist_ok=True)

    for col in adata.obs.columns: 
        adata.obs[col] = adata.obs[col].astype("category")

    adata = add_gene_binary_columns(adata, genes_list)

    # --- FIX 1: ensure ALL needed color columns exist ---
    for col in cols_to_plot:
        adata, new_gene_colors, new_cell_colors = set_obs_colors(adata, palette, cell_type_colors, gene_color_map)

    gene_color_map = {**gene_color_map, **new_gene_colors}
    cell_type_colors = {**cell_type_colors, **new_cell_colors}

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
        
        if plot_always == False and os.path.exists(umap_out_path):
            continue

        legend_loc = "on data" if col in ["leiden"] else "right margin"
        title = f"{sample_id}\n{col}"

        # --- FIX 3: gene-specific coloring ---
        palette_arg = None

        if col.endswith("_binary"):
            gene = col.replace("_binary", "")

            if gene in gene_counts:
                n_pos, n_total = gene_counts[gene]
                title = f"{sample_id}\n{gene} ({n_pos}/{n_total} cells positive)"

            # Apply gene-specific color if available
            if gene in gene_color_map:
                palette_arg = ["lightgrey", gene_color_map[gene]]

        n_categories = adata.obs[col].nunique()

        width = 5 + min(n_categories * 0.4, 8)
        fig, ax = plt.subplots(
            figsize=(width, 5) if legend_loc == "right margin" else (5, 5)
        )

        sc.pl.umap(
            adata,
            color=col,
            legend_loc=legend_loc,
            frameon=False,
            ax=ax,
            title=title,
            show=False,
            palette=palette_arg  # <-- key fix
        )

        if legend_loc == "right margin":
            plt.subplots_adjust(right=0.75)

        plt.tight_layout()
        plt.savefig(umap_out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return adata, gene_color_map, cell_type_colors
            

# One driver function (what you actually run)
import os
import scanpy as sc
import numpy as np
import pandas as pd
import json


def show_cta_genes(project_name, adata, figures_dir, tables_dir,cta_genes_file_path, groupby, cta_genes_fixed=None):

    cta_tables_save_dir = os.path.join(tables_dir, "cta_analysis")
    os.makedirs(cta_tables_save_dir, exist_ok=True)
    cta_figures_save_dir = os.path.join(figures_dir, "cta_analysis")
    os.makedirs(cta_figures_save_dir, exist_ok=True)
    sc.settings.figdir = cta_figures_save_dir

    if cta_genes_fixed is None:
        cta_genes_df = pd.read_csv(cta_genes_file_path)
        cta_genes_all = cta_genes_df["Family member"].unique().tolist()
    else:
        cta_genes_all = cta_genes_fixed
    
    cta_genes_in = []
    gene_stats = []

    for g in cta_genes_all:
        if g not in adata.raw.var_names:
            continue
        else:
            cta_genes_in.append(g)

        vals = adata.raw[:, g].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray().flatten()
        else:
            vals = vals.flatten()

        frac = np.mean(vals > 0)
        mean_expr = np.mean(vals)

        gene_stats.append({"gene": g, "fraction_expressing": frac, "mean_expression": mean_expr})

    df_stats = pd.DataFrame(gene_stats)
    df_stats.to_csv(os.path.join(cta_tables_save_dir, f"{project_name}_cta_gene_stats.csv"), index=False)

    df_expr = sc.get.obs_df(adata, keys=[groupby] + cta_genes_in, use_raw=True)
    df_mean = df_expr.groupby(groupby).mean()
    df_mean.to_csv(os.path.join(cta_tables_save_dir, f"{project_name}_cta_mean_expression.csv"))

    df_frac = df_expr.copy()
    for g in cta_genes_in:
        df_frac[g] = df_frac[g] > 0

    df_frac = df_frac.groupby(groupby).mean()
    df_frac.to_csv(os.path.join(cta_tables_save_dir, f"{project_name}_cta_fraction_expression.csv"))

    df_long = df_expr.melt(id_vars=groupby, var_name="gene", value_name="expression")
    df_long.to_csv(os.path.join(cta_tables_save_dir, f"{project_name}_cta_long_format.csv"), index=False)

    
    cta_genes_plot  = cta_genes_in
    #cta_genes_plot  = [g for g in cta_genes_in if g in adata.var_names]
   
    print("Final CTA genes:", cta_genes_plot)

    sc.pl.dotplot(adata, var_names=cta_genes_plot, groupby=groupby, cmap="Purples", use_raw=True, standard_scale=None, title=f"{project_name} CTA genes", save=f"{project_name}_cta_dotplot.png")
    sc.pl.umap(adata, color=[groupby], frameon=False, save=f"{project_name}_{groupby}_umap.png")
    sc.pl.umap(adata, color= cta_genes_plot, use_raw=True, ncols=3, cmap="Purples", vmax="p99", frameon=False, save=f"{project_name}_cta_umap.png")

    return cta_genes_plot
