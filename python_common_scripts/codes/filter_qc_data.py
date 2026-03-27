# python_common_scripts/codes/filter_qc_data.py

#### Libraries ###
import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc

import anndata as ad
import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.sparse import issparse
from scipy.io import mmread

sc.settings.autosave = True  # Saves plots to figdir automatically
sc.settings.autoshow = False # Prevents plots from displaying inline

import matplotlib.pyplot as plt

#### Custom ####

from config import *
#from codes.cell_annotation import *
#from codes.module_scoring import cell_cycle_scoring, program_scoring

#### Functions ####

def filter_data(sample_id, adata_init, adata_filtered_path, important_genes, cell_cycle_genes_file_path, qc_save_dir, full_filter = 0.05, relaxed_filter = 0.02):
    
    #  Setup 
    violin_dir = os.path.join(qc_save_dir, "violin")
    scatter_dir = os.path.join(qc_save_dir, "scatter")

    os.makedirs(qc_save_dir, exist_ok=True)
    os.makedirs(violin_dir, exist_ok=True)
    os.makedirs(scatter_dir, exist_ok=True)

    full_filter = float(full_filter)
    relaxed_filter = float(relaxed_filter)
    adata = adata_init.copy()

    # Remove cells with no genes and genes with no expression
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)

    # Calculate and Plot QC Metrics
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB(?!P)", regex=True)

    # Calculate the QC metrics for the defined sets
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=False)

    # Plot the QC metrics BEFORE applying cell filters
    qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    
    sc.settings.figdir = violin_dir
    sc.pl.violin(adata, qc_cols, save=f"_{sample_id}_QC.png", jitter=0.4, multi_panel=True)
    
    sc.settings.figdir = scatter_dir
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save=f"_{sample_id}_QC.png")

    # Create a boolean mask to identify and remove low-quality cells

    if adata.n_obs > 50:

        cell_mask = (
            (adata.obs["n_genes_by_counts"] >= adata.obs["n_genes_by_counts"].quantile(full_filter)) &
            (adata.obs["total_counts"] >= adata.obs["total_counts"].quantile(full_filter)) &
            (adata.obs["pct_counts_mt"] <= adata.obs["pct_counts_mt"].quantile(1 - full_filter))
        )
        adata = adata[cell_mask].copy()

    else:

        cell_mask = (
            (adata.obs["n_genes_by_counts"] >= adata.obs["n_genes_by_counts"].quantile(relaxed_filter)) &
            (adata.obs["total_counts"] >= adata.obs["total_counts"].quantile(relaxed_filter)) &
            (adata.obs["pct_counts_mt"] <= adata.obs["pct_counts_mt"].quantile(1 - relaxed_filter))
        )
        adata = adata[cell_mask].copy()

    

    # Create a mask for genes expressed in at least 3 cells
    genes_expressed_mask = sc.pp.filter_genes(adata, min_cells=3, inplace=False)[0]

    # Create a mask for the important genes
    important_set = {g.upper() for g in important_genes}
    important_genes_mask = adata.var_names.str.upper().isin(important_set)

    final_gene_mask = genes_expressed_mask | important_genes_mask

    # Apply this final, combined mask to the data just once
    adata = adata[:, final_gene_mask].copy()

    cell_cycle_scoring(sample_id, adata, cell_cycle_genes_file_path, qc_save_dir)

    # Store the raw counts in a layer before any normalization
    adata.layers["counts"] = adata.X.copy()
    adata.write(adata_filtered_path)

    return adata

def cell_cycle_scoring(sample_id, adata, cell_cycle_genes_file_path, qc_save_dir):

    cell_cycle_dir = os.path.join(qc_save_dir, "cell_cycle")
    os.makedirs(cell_cycle_dir, exist_ok=True)

    cell_cycle_genes_df = pd.read_csv(cell_cycle_genes_file_path, sep="\t")

    s_genes = cell_cycle_genes_df.loc[cell_cycle_genes_df.phase == "G1-S", "gene"].tolist()
    s_genes = [g for g in s_genes if g in adata.var_names]

    g2m_genes = cell_cycle_genes_df.loc[cell_cycle_genes_df.phase == "G2-M", "gene"].tolist()
    g2m_genes = [g for g in g2m_genes if g in adata.var_names]

    print(f"Found {len(s_genes)} G1-S genes and {len(g2m_genes)} G2-M genes in {sample_id}")
    if len(s_genes) == 0 or len(g2m_genes) == 0:
        print("Insufficient matching cell cycle genes found.")
        return 
    
    else:
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes,  g2m_genes=g2m_genes)

        sc.settings.figdir =cell_cycle_dir
        sc.pl.violin( adata,["S_score", "G2M_score"], groupby="phase", jitter=0.4, multi_panel=True,show=False, save=f"_{sample_id}_cell_cycle_.png")
        
def run_qc(adata_filtered, adata_qc_path, important_genes):
    adata = adata_filtered.copy()
    
    # always restore raw counts into X
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        
    
    # These two steps correctly prepare your data in adata.X
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Save the log-normalized data before filtering to HVGs
    adata.raw = adata.copy()
    
    # Find HVGs using the log-normalized data in adata.X by removing the 'layer' argument.
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
    
    # Force important genes to be included as "highly variable"
    important_set = {g.upper() for g in important_genes}
    gene_series = adata.var_names.astype(str)
    adata.var.loc[gene_series.str.upper().isin(important_set),"highly_variable"] = True
    
    # Filter the data to only the highly variable genes
    adata = adata[:, adata.var["highly_variable"]].copy()
    
    # Continue with downstream analysis
    sc.pp.scale(adata, max_value=10)
    max_pcs = max(2, min(adata.n_obs, adata.n_vars)-1)
    n_pcs = min(30, max_pcs)
    sc.tl.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0)

    adata.write(adata_qc_path)
    return adata
        