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


def get_ranked_genes(adata, sample_id, tables_dir, figures_dir, always_rank, marker_genes_dict, save_name=None):

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


def plot_marker_genes_leiden(adata,adata_path, cell_type_colors, gene_color_map,leiden_res=0.50):

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
    adata.write(adata_path)
    return adata

#adata_male_path = os.path.join(working_adata_dir, f"{project_name}_male_processed.h5ad")
#adata_female_path = os.path.join(working_adata_dir, f"{project_name}_female_processed.h5ad")

import json


def rank_genes_for_publication(adata_female_path, project_name, final_leiden_res, cluster_celltype_mapping_json, cta_genes_file_path, figures_dir):

    annotation_figures_dir = os.path.join(figures_dir, "annotation")
    os.makedirs(annotation_figures_dir, exist_ok=True)
    
    sc.settings.autosave = False
    sc.settings.autoshow = False
    sc.settings.figdir = annotation_figures_dir
    sc.settings.set_figure_params(dpi=300)
    
    final_leiden_res = str(final_leiden_res)
    final_key = f"leiden_res_split_{final_leiden_res}"
    with open(cluster_celltype_mapping_json, "r") as f:
        data = json.load(f)

    # load both resolutions
    map_3 = data[f"{project_name}_female"]["leiden_res_3.0"]["cluster_to_celltype"]
    map_2 = data[f"{project_name}_female"]["leiden_res_2.0"]["cluster_to_celltype"]
    
    cta_genes_umap_save_path = f"_{project_name}_female_cta_genes_{final_leiden_res}.png"
    cta_genes_dotplot_save_path = f"{project_name}_female_cta_genes_{final_leiden_res}.png"
    leiden_annot_genes_dotplot_save_path = f"{project_name}_female_annot_leiden_{final_leiden_res}.png"
    ranked_genes_dotplot_save_path = f"{project_name}_female_ranked_genes_leiden_{final_leiden_res}.png"
    umap_out_path = f"_{project_name}_female_ranked_genes_leiden_{final_leiden_res}.png"

    adata = sc.read_h5ad(adata_female_path)
    
    sc.pp.neighbors(adata, use_rep="X_pca", n_pcs=30)
    sc.tl.umap(adata)
    # Clustering
    for res in [0.5, 1.0, 2.0,3.0]:
        sc.tl.leiden(adata, key_added=f"leiden_res_split_{res:3.1f}", resolution=res, flavor="igraph")
    
    adata.write(adata_female_path)
    
    # 1. Get markers per  cluster
    sc.tl.rank_genes_groups(adata, groupby=final_key,method="wilcoxon", use_raw=True)

    #2. Extract top markers
    markers = sc.get.rank_genes_groups_df(adata, group=None)

    # top 10 per cluster
    top_markers = (markers.groupby("group", observed=False).head(10).reset_index(drop=True))
    top_markers = top_markers[top_markers["names"].notna()]
    top_markers_leiden_dict = top_markers.groupby("group", observed=False)["names"].apply(lambda x: x.tolist()[:5]).to_dict()
    
    sc.tl.dendrogram(adata, groupby=final_key)

    # 5. Dotplot
    sc.pl.dotplot(adata,  top_markers_leiden_dict, groupby=final_key, dendrogram=True, standard_scale="var", cmap="Purples", var_group_rotation=0, dot_max=0.6, dot_min=0.1, figsize=(15, 6), show=False, save=ranked_genes_dotplot_save_path)
    
    #6. Map to cells
   
    adata.obs["celltype_2"] = adata.obs["leiden_res_split_2.0"].map(map_2) #map 2.0 labels to all cells
    adata.obs["celltype_3"] = adata.obs["leiden_res_split_3.0"].map(map_3) #map 3.0 labels (partial, more specific)
    #Use 3.0 only where it’s defined, otherwise keep 2.0:

    adata.obs["final_cell_type"] = adata.obs["celltype_3"].where( adata.obs["celltype_3"].notna() & (adata.obs["celltype_3"] != ""), adata.obs["celltype_2"]).fillna("Unknown")
    final_cluster_celltype_dict = {map_3.get(k, k): v for k, v in top_markers_leiden_dict.items()}

    adata.obs[f"final_cell_type_{final_leiden_res}"] = adata.obs["final_cell_type"]
    
    sc.pl.dotplot(adata, final_cluster_celltype_dict, groupby=final_key, dendrogram=True, standard_scale="var", cmap="Purples", var_group_rotation=0, dot_max=0.6, dot_min=0.1, figsize=(15, 6), show=False, save=leiden_annot_genes_dotplot_save_path)
    
    # 7. Visualize on umap
    sc.pl.umap(adata, color = [f"final_cell_type_{final_leiden_res}",final_key], title = ["Female Cell Type Annotation", f"Leiden Clustering (resolution {final_leiden_res})"], legend_loc = "on data", frameon = False, save = umap_out_path)

    # 8. Visualize cta genes
    cta_genes_df = pd.read_csv(cta_genes_file_path)
    cta_genes = cta_genes_df["Family member"].to_list()
    cta_genes_for_plot = [g for g in cta_genes if g in adata.var_names]

    sc.pl.dotplot(adata, cta_genes_for_plot, groupby=f"final_cell_type_{final_leiden_res}", dendrogram=True, standard_scale="var", cmap="Purples", var_group_rotation=0, dot_max=0.6,dot_min=0.1, figsize=(15, 6),show=False, save=cta_genes_dotplot_save_path)
    sc.pl.umap(adata, color = [f"final_cell_type_{final_leiden_res}"] + cta_genes_for_plot, title= ["Female Cell Type Annotation"], cmap = "Purples", frameon = False, ncols=4, legend_loc = "on data", save = cta_genes_umap_save_path)

    adata.write(adata_female_path)
    
    return adata

def markers_leiden_embryos_mixed(adata, subproject, top_markers_path, final_leiden_res, figures_dir):
    deg_dir = os.path.join(figures_dir, "deg_analysis")
    os.makedirs(deg_dir, exist_ok=True)
    sc.settings.figdir = deg_dir

    leiden_split = f"leiden_res_split_{final_leiden_res}"
    leiden_celltype_split = f"final_celltype_{final_leiden_res}"
    leiden_cellstate_split = f"final_cellstate_{final_leiden_res}"

    adata = adata.copy()

    if os.path.exists(top_markers_path):
        top_markers = pd.read_csv(top_markers_path)
        top5 = top_markers.groupby("group")["names"].head(5).tolist()
    else:
        sc.tl.rank_genes_groups(adata, groupby=leiden_split, method="wilcoxon", use_raw=True)

        markers = sc.get.rank_genes_groups_df(adata, group=None)
        top_markers = markers.groupby("group").head(10).reset_index(drop=True)
        top_markers.to_csv(top_markers_path)

        top5 = top_markers.groupby("group")["names"].head(5).tolist()

        sc.tl.dendrogram(adata, var_names=top5, groupby=leiden_split)
        sc.pl.dendrogram(adata, groupby=leiden_split, save=f"_{subproject}_{leiden_split}.png")
        sc.pl.dotplot(adata, var_names=top5, groupby=leiden_split, standard_scale="var", cmap="Purples", dendrogram=True, save=f"_{subproject}_{leiden_split}_top5.png")

    sc.pl.umap(
        adata,
        color=[leiden_split, leiden_celltype_split, leiden_cellstate_split],
        legend_loc="on data",
        frameon=False,
        save=f"_{subproject}_{leiden_split}.png"
    )