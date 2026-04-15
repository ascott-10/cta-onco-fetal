# python_common_scripts/codes/differential_gene_expression.py

#### Libraries ###
import os
import glob
import json
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


def get_ranked_genes(adata, adata_ranked_path, subproject_name, json_annotations_path, leiden_res_list, tables_dir, figures_dir, always_rank=True):

    markers_save_dir = os.path.join(figures_dir, "deg_analysis")
    os.makedirs(markers_save_dir, exist_ok=True)
    sc.settings.figdir = markers_save_dir

    if os.path.exists(json_annotations_path):
        with open(json_annotations_path) as f:
            annotations = json.load(f)
        sub_ann = annotations.get(subproject_name, {})
    else:
        print("Annotation file not found")
        sub_ann = {}

    print("debug sub_ann keys:", sub_ann.keys())
    leiden_cols = []
    celltype_cols = []
    leiden_cols_titles = []

    for res in leiden_res_list:

        key_string = f"leiden_res_{res}"
        celltype_col = f"celltype_{key_string}"

        if key_string not in adata.obs:
            continue

        print(f"debug {key_string} in JSON:", key_string in sub_ann)
        
        adata.obs[key_string] = adata.obs[key_string].astype(str).astype("category")

        top_markers_path = os.path.join(tables_dir, f"{subproject_name}_{key_string}.csv")

        if os.path.exists(top_markers_path):
            top_markers = pd.read_csv(top_markers_path)
        else:
            top_markers = None

        use_raw = adata.raw is not None

        if not (os.path.exists(top_markers_path) and not always_rank):
            sc.tl.rank_genes_groups(adata, groupby=key_string, method="wilcoxon", use_raw=use_raw)
            markers = sc.get.rank_genes_groups_df(adata, group=None)
            top_markers = markers.groupby("group").head(10).reset_index(drop=True)
            top_markers.to_csv(top_markers_path, index=False)

        top_markers_dict = top_markers.groupby("group")["names"].apply(list).to_dict()

        top5 = top_markers.groupby("group")["names"].apply(lambda x: x.head(5).tolist()).to_dict()

        sc.tl.dendrogram(adata, var_names=[g for v in top5.values() for g in v], groupby=key_string)
        sc.pl.dendrogram(adata, groupby=key_string)
        sc.pl.dotplot(adata, var_names=top5, groupby=key_string, standard_scale="var", cmap="Purples", dendrogram=True, title=f"{subproject_name} \n Top Markers ({key_string})", save=f"{subproject_name}_top_markers_{key_string}.png")

        print(f"\n=== {key_string} ===")
        for k, v in top_markers_dict.items():
            print(k, v)

        leiden_ann = sub_ann.get(key_string, {})
        print("debug leiden_ann:", leiden_ann)

        cluster_to_celltype_dict = leiden_ann.get("cluster_to_celltype", {})
        cluster_to_celltype_dict = {str(k): v for k, v in cluster_to_celltype_dict.items()} #force to string

        print("debug cluster_to_celltype_dict:", cluster_to_celltype_dict)

        print("debug adata clusters:", adata.obs[key_string].unique())
        print("debug mapping keys:", cluster_to_celltype_dict.keys())

        if cluster_to_celltype_dict:

            print(f"debug Creating column: {celltype_col}")
            adata.obs[celltype_col] = adata.obs[key_string].astype(str).map(cluster_to_celltype_dict).fillna("Unknown").astype("category")
            mapped = adata.obs[key_string].astype(str).map(cluster_to_celltype_dict)
            print("debug", adata.obs[celltype_col].value_counts())

            print("debug mapped unique:", mapped.unique())
            print("debug num unmapped:", mapped.isna().sum())


            sc.tl.dendrogram(adata, var_names=[g for v in top5.values() for g in v], groupby=celltype_col)
            sc.pl.dotplot(adata, var_names=top5, groupby=celltype_col, standard_scale="var", cmap="Purples", dendrogram=True, title=f"{subproject_name} \n Top Markers (celltype)", save=f"{subproject_name}_top_markers_celltype_{key_string}.png")

        celltype_cols.append(celltype_col)
        leiden_cols.append(key_string)
        leiden_cols_titles.append(f"Leiden Clustering \n Resolution {res}")

    adata.write(adata_ranked_path)

    return adata, leiden_cols, leiden_cols_titles, celltype_cols

