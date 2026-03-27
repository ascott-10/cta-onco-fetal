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
            

def cta_genes_expression(sample_id, adata, tables_dir, figures_dir, cta_genes_file_path, top_n, plot_always=False):
    """No return"""

    cta_table_save_dir = os.path.join(tables_dir, "cta_genes")
    os.makedirs(cta_table_save_dir, exist_ok = True)
    cta_figure_save_dir = os.path.join(figures_dir, "cta_genes")
    os.makedirs(cta_figure_save_dir, exist_ok = True)

    all_expression_save_path = os.path.join(cta_table_save_dir, f"{sample_id}_cta_genes_all_expression.csv")
    dotplot_save_path = os.path.join(cta_figure_save_dir,f"{sample_id}_cta_genes_dotplot.png")
    
    if os.path.exists(all_expression_save_path):
        df_save = pd.read_csv(all_expression_save_path)
    else:
    
        # Read in CTA genes list
        df = pd.read_csv(cta_genes_file_path)
        df.columns = ["gene_group", "gene"]

        all_genes_list = sorted(df["gene"].unique())
        
        # 1. Make table for each gene, celltype, mean expression, fraction expression

        

        rows = []

        for gene in all_genes_list:

            if gene not in adata.var_names:
                continue

            expr = adata[:, gene].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()

            for cell_type in adata.obs["predicted_cell_type"].unique():
            
                mask = adata.obs["predicted_cell_type"] == cell_type
                expr_ct = expr[mask.values]

                n_ct = mask.sum()
                n_pos = (expr_ct > 0).sum()
                mean_expr = expr_ct.mean() if n_ct > 0 else 0.0
                frac_pos = n_pos / n_ct if n_ct > 0 else 0.0

                rows.append({
                    "sample_id": sample_id,
                    "predicted_cell_type": cell_type,
                    "gene": gene,
                    "n_cells": int(n_ct),
                    "n_positive_cells": int(n_pos),
                    "frac_positive_cells": float(frac_pos),
                    "mean_expression": float(mean_expr),
                })

        df_save = pd.DataFrame(rows)
        df_save.to_csv(all_expression_save_path, index=False)
    
    if (not os.path.exists(dotplot_save_path) or plot_always == True):
       
        # 2. To visualize CTA genes across dataset
        gene_dataset = (df_save.groupby(["gene", "predicted_cell_type"])
                        .agg({"mean_expression": "mean", "frac_positive_cells": "mean"})
                        .reset_index())
        
        
        gene_order = (gene_dataset.groupby("gene")["mean_expression"].max().sort_values(ascending=False).index) # Order by max expression

        # Create df for dotplot
        dotplot_df = gene_dataset.copy()
        
        gene_order_top = gene_order[:top_n]

        dotplot_df = dotplot_df[dotplot_df["gene"].isin(gene_order_top)]
        dotplot_df["gene"] = pd.Categorical(
            dotplot_df["gene"],
            categories=gene_order_top,
            ordered=True
        )

        # Cluster for matrix - mean expression
        cluster_df = (dotplot_df.pivot(index="predicted_cell_type",columns="gene",values="mean_expression").fillna(0))

        # Cluster cell types 
        X = cluster_df.values
        Z = linkage(X, method="average", metric="correlation")
        cell_order = cluster_df.index[leaves_list(Z)]

        # Apply order
        dotplot_df["predicted_cell_type"] = pd.Categorical(dotplot_df["predicted_cell_type"],categories=cell_order,ordered=True)

        # Make the Dotplot 
        plt.figure(figsize=(top_n * 0.3, len(cell_order) * 0.5))

        sns.scatterplot(
            data=dotplot_df,
            x="gene",
            y="predicted_cell_type",
            size="frac_positive_cells",
            hue="mean_expression",
            sizes=(10, 200),
            palette="Purples",
            hue_norm=(0, dotplot_df["mean_expression"].quantile(0.95)),
            edgecolor="none"
        )

        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.yticks(fontsize=9)

        plt.xlabel("")
        plt.ylabel("")
        plt.title(f"{sample_id}\nCTA Gene Expression", fontsize=12)

        plt.grid(True, axis="x", linewidth=0.3, alpha=0.3)
        plt.grid(False, axis="y")

        plt.legend([], [], frameon=False)

        plt.tight_layout()

        plt.savefig(dotplot_save_path, dpi=300, bbox_inches="tight")
        plt.close()
    
    return df_save

    
def cta_genes_expression_all_samples(project_name):

    df = pd.read_csv(f"all_results/tables/cta_genes/{project_name}_cta_gene_expression.csv")

    dotplot_save_dir = "all_results/figures/cta_genes"
    os.makedirs(dotplot_save_dir, exist_ok = True)
    df = pd.read_csv(f"all_results/tables/cta_genes/{project_name}_cta_gene_expression.csv")
    # Collapse to sample x gene
    sample_gene_df = (
        df.groupby(["sample_id", "gene"])
        .apply(lambda x: pd.Series({
            "n_cells": x["n_cells"].sum(),
            "n_positive_cells": x["n_positive_cells"].sum(),
            "frac_positive_cells": x["n_positive_cells"].sum() / x["n_cells"].sum(),
            "mean_expression": (x["mean_expression"] * x["n_cells"]).sum() / x["n_cells"].sum()
        }))
        .reset_index()
    )

    # Pick genes
    top_n = 25

    gene_order = (
        sample_gene_df.groupby("gene")["mean_expression"]
        .max()
        .sort_values(ascending=False)
        .index
    )

    gene_order_top = gene_order[:top_n]

    plot_df = sample_gene_df[sample_gene_df["gene"].isin(gene_order_top)].copy()

    # Set ordering
    plot_df["gene"] = pd.Categorical(
        plot_df["gene"],
        categories=gene_order_top,
        ordered=True
    )

    # Cluster samples - optional but useful
    from scipy.cluster.hierarchy import linkage, leaves_list

    cluster_df = (
        plot_df.pivot(index="sample_id", columns="gene", values="mean_expression")
        .fillna(0)
    )

    Z = linkage(cluster_df.values, method="average", metric="correlation")
    sample_order = cluster_df.index[leaves_list(Z)]

    plot_df["sample_id"] = pd.Categorical(
        plot_df["sample_id"],
        categories=sample_order,
        ordered=True
    )

    import matplotlib.pyplot as plt
    import seaborn as sns

    height = max(4, len(cell_order) * 0.6)

    plt.figure(figsize=(top_n * 0.35, height))

    sns.scatterplot(
        data=plot_df,
        x="gene",
        y="sample_id",
        size="frac_positive_cells",
        hue="mean_expression",
        sizes=(10, 200),
        palette="Purples",
        edgecolor="none",
        hue_norm=(0, plot_df["mean_expression"].quantile(0.95))
    )

    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=9)

    plt.xlabel("")
    plt.ylabel("")
    plt.title(f"{project_name}\nCTA Gene Expression Across Samples")

    plt.legend([], [], frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(dotplot_save_dir, f"{project_name}_cta_all_samples_dotplot.png"),dpi=300, bbox_inches="tight")
    plt.close()