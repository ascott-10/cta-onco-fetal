# python_common_scripts/codes/differential_gene_expression.py

#### Libraries ###
import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc

sc.settings.autosave = True  # Saves plots to figdir automatically
sc.settings.autoshow = False # Prevents plots from displaying inline

# Needed for some plotting
import matplotlib.pyplot as plt

import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread

import matplotlib.pyplot as plt

#### Custom ####

from config import *
from codes.cell_annotation import *


def get_ranked_genes(adata, sample_id, tables_dir, figures_dir, always_rank=True, marker_genes_dict=CUSTOM_MARKER_GENES_DICT, save_name=None):

    deg_tables_dir = os.path.join(tables_dir, "deg_analysis")
    deg_figures_dir = os.path.join(figures_dir, "deg_analysis")
    os.makedirs(deg_figures_dir, exist_ok=True)
    os.makedirs(deg_tables_dir, exist_ok=True)

    if save_name == None:
        ranked_genes_dotplot_save_path = f"_{sample_id}_ranked_genes_2_dotplot.pdf"
        cluster_ranked_genes_df_save_path = os.path.join(deg_tables_dir, f"{sample_id}_top_genes_per_leiden_cluster.csv")
    else:
        ranked_genes_dotplot_save_path = f"_{sample_id}_{save_name}_ranked_genes_2_dotplot.pdf"
        cluster_ranked_genes_df_save_path = os.path.join(deg_tables_dir, f"{sample_id}_{save_name}_top_genes_per_leiden_cluster.csv")

    # --- work on copy (fix)
    ad = adata.copy()

    # Only use highly variable genes
    ad.X = np.nan_to_num(ad.X, nan=0.0, posinf=0.0, neginf=0.0)

    if "highly_variable" not in ad.var.columns:
        sc.pp.highly_variable_genes(
            ad,
            layer="counts",
            n_top_genes=2000,
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            flavor="seurat_v3"
        )

    ad = ad[:, ad.var["highly_variable"]].copy()

    # 1. Create df of top genes in each cluster
    

    if os.path.exists(cluster_ranked_genes_df_save_path) and not always_rank:
        cluster_markers_df = pd.read_csv(cluster_ranked_genes_df_save_path)

    else:
        sc.tl.rank_genes_groups(ad, "leiden", mask_var="highly_variable", method="t-test", use_raw=False)

        cluster_dfs = []

        for group in sorted(ad.obs["leiden"].unique()):
            df = sc.get.rank_genes_groups_df(ad, group)

            df_filtered = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 1)]
            if df_filtered.shape[0] >= 5:
                df = df_filtered

            df = df[["names", "logfoldchanges", "pvals_adj"]].head(10).copy()
            df["leiden_cluster"] = group
            df["rank"] = range(1, len(df) + 1)

            cluster_dfs.append(df)

        cluster_markers_df = pd.concat(cluster_dfs, ignore_index=True)
        cluster_markers_df.to_csv(cluster_ranked_genes_df_save_path, index=False)

    # 2. Dotplot
    

    marker_genes_filtered = {
        k: [g for g in v if g in ad.var_names]
        for k, v in marker_genes_dict.items()
    }

    marker_genes_filtered = {
        k: v for k, v in marker_genes_filtered.items() if len(v) > 0
    }

    clusters = list(ad.obs["leiden"].unique())

    sc.tl.dendrogram(ad, groupby="leiden")

    sc.settings.figdir = deg_figures_dir

    sc.pl.dotplot(
        ad,
        marker_genes_filtered,
        groupby="leiden",
        categories_order=clusters,
        dendrogram=True,
        standard_scale="var",
        cmap="Purples",
        var_group_rotation=0,
        dot_max=0.6,
        dot_min=0.1,
        figsize=(15, 6),
        show=False,
        save=ranked_genes_dotplot_save_path  # <-- fixed
    )
    
    return adata


def plot_marker_genes(adata, dataset_name, marker_dict, figures_dir, new_col, leiden_res=0.50):

    umap_dir = os.path.join(figures_dir, "marker_genes", "umap")
    os.makedirs(umap_dir, exist_ok=True)
    umap_save_path = os.path.join(umap_dir, f"{dataset_name}_new_annotated_umap.png")

    dotplot_dir = os.path.join(figures_dir, "marker_genes", "dotplot")
    os.makedirs(dotplot_dir, exist_ok=True)
    dotplot_save_path = os.path.join(dotplot_dir, f"{dataset_name}_new_annotated_dotplot.png")

    # Compute neighbors + UMAP if not already present
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata)
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    # Leiden
    leiden_key = f"leiden_res_{leiden_res:4.2f}"
    sc.tl.leiden(adata, key_added=leiden_key, resolution=leiden_res, flavor="igraph")

    

    # Ensure categorical
    color_list = [leiden_key, "sample_id", "predicted_cell_type"]
    for col in color_list:
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype("category")

    adata, gene_color_map, cell_type_colors = set_obs_colors(
        adata,
        palette="glasbey",
        cell_type_colors=cell_type_colors,
        gene_color_map=gene_color_map
    )

    # --- UMAP ---
    sc.pl.umap(
        adata,
        color=color_list,
        ncols=len(color_list),
        frameon=False,
        show=False
    )

    fig = plt.gcf()
    axes = fig.axes

    for ax, col in zip(axes, color_list):
        if col == leiden_key:
            leg = ax.get_legend()
            if leg:
                leg.set_bbox_to_anchor((0.5, 0.5))
        ax.set_title(f"{dataset_name}\n{col}")

    plt.tight_layout()
    fig.savefig(umap_save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def stop():

    # Apply mapping
    adata, cluster_order = apply_celltype_mapping(
        adata,
        mapping_file=CELL_TYPE_MAPPINGS,
        dataset_name=dataset_name,
        new_col=new_col,
        leiden_key=leiden_key
    )

    
    # --- DOTPLOT ---
    dot = sc.pl.dotplot(
        adata,
        var_names=marker_dict,
        groupby=leiden_key,
        categories_order=cluster_order,
        standard_scale="var",
        cmap="Purples",
        return_fig=True
    )

    dot.make_figure()
    ax = dot["mainplot_ax"]

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center")

    dot["fig"].suptitle(
        f"{dataset_name}\nMarker Genes for Key Cell Types",
        y=1.02
    )

    dot["fig"].tight_layout()
    dot["fig"].subplots_adjust(top=0.88)

    dot["fig"].savefig(dotplot_save_path, dpi=300, bbox_inches="tight")
    plt.close(dot["fig"])

    return adata