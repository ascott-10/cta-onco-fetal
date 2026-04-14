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


def filter_data_embryos_mixed(project_name, adata_init, adata_filt_path, important_genes, qc_save_dir):
    """For consistency:
        counts" = true raw data
        adata.X = log-transformed
        adata.raw = clean log reference
        """
    # Set up dirs
    violin_dir = os.path.join(qc_save_dir, "violin")
    
    os.makedirs(qc_save_dir, exist_ok=True)
    os.makedirs(violin_dir, exist_ok=True)

    adata = adata_init.copy()

    adata.obs["nonzero_counts"] = np.array((adata.X > 0).sum(axis=1)).flatten()  # Because original data was CPM or TPM:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-") # mitochondrial genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL")) # ribosomal genes
    adata.var["hb"] = adata.var_names.str.startswith("HB") & ~adata.var_names.str.startswith("HBP") # hemoglobin genes

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True)

    qc_cols = ["n_genes_by_counts", "nonzero_counts", "pct_counts_mt"] # Because original data was CPM or TPM, normully use total counts not nonzero counts

    sc.settings.figdir = qc_save_dir
    sc.pl.violin(adata, qc_cols, save=f"_{project_name}_QC.png", jitter=0.4, multi_panel=True)
 
    # Filter for cells and genes here
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)

    # Preserve raw counts BEFORE log transform
    if "raw_counts" in adata.layers:
        adata.layers["counts"] = adata.layers["raw_counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    # Create a mask for genes expressed in at least 3 cells
    genes_expressed_mask = sc.pp.filter_genes(adata, min_cells=3, inplace=False)[0]

    # Create a mask for the important genes
    important_set = {g.upper() for g in important_genes}
    important_genes_mask = adata.var_names.str.upper().isin(important_set)

    final_gene_mask = genes_expressed_mask | important_genes_mask

    # Apply this final, combined mask to the data just once
    adata = adata[:, final_gene_mask].copy()

    

    adata.write(adata_filt_path)

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

def run_qc_embryos_mixed(project_name, adata_filtered, adata_qc_path, important_genes, qc_save_dir, cell_cycle_genes_file_path):
    """For consistency:
        counts" = true raw data
        adata.X = log-transformed
        adata.raw = clean log reference
        """
    pca_dir = os.path.join(qc_save_dir, "pca")
    os.makedirs(pca_dir, exist_ok=True)
    sc.settings.figdir = pca_dir
    
    n_before = adata_filtered.n_obs
    adata = adata_filtered.copy()


     # Remove doublets FIRST
    if "predicted_doublet" in adata.obs:
        adata = adata[~adata.obs["predicted_doublet"]].copy()

    n_after = adata.n_obs
    print(f"Removed {n_before - n_after} doublets in project {project_name}")


    # Restore raw counts into X
    if "raw_counts" in adata.layers:
        adata.X = adata.layers["raw_counts"].copy()
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    # Log transform
    sc.pp.log1p(adata)

    # Store log-normalized version
    adata.raw = adata.copy()

    cell_cycle_scoring(project_name, adata, cell_cycle_genes_file_path, qc_save_dir)

    
    # HVGs
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample_id", flavor="seurat")
    sc.pl.highly_variable_genes(adata, save=f"_{project_name}_HVGs.png")
    
    # Force important genes into HVGs
    important_set = {g.upper() for g in important_genes}
    gene_series = adata.var_names.astype(str)
    adata.var.loc[gene_series.str.upper().isin(important_set), "highly_variable"] = True
    
    # Subset to HVGs
    adata = adata[:, adata.var["highly_variable"]].copy()

    # ADD THIS BLOCK
    adata.obs['total_counts'] = pd.to_numeric(adata.obs['total_counts'], errors='coerce')

    if ("S_score" in adata.obs) and ("G2M_score" in adata.obs):
        sc.pp.regress_out(adata, ['total_counts', 'S_score', 'G2M_score'])
    else:
        sc.pp.regress_out(adata, ['total_counts'])

    # PCA
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save=f"_{project_name}_PCA_variance_ratio.png")

    # Neighbors
    sc.pp.neighbors(adata, n_pcs=30)

    # Clustering
    for res in [0.5, 1.0, 2.0, 3.0]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:3.1f}", resolution=res, flavor="igraph")

    # UMAP
    sc.tl.umap(adata)

    adata.write(adata_qc_path)
    
    return adata

def filter_data_ovarian(subproject_name, adata_init, adata_filt_path, important_genes, qc_save_dir, full_filter=0.05, relaxed_filter=0.02):

    violin_dir = os.path.join(qc_save_dir, "violin")
    scatter_dir = os.path.join(qc_save_dir, "scatter")

    os.makedirs(qc_save_dir, exist_ok=True)
    os.makedirs(violin_dir, exist_ok=True)
    os.makedirs(scatter_dir, exist_ok=True)

    adata = adata_init.copy()

    full_filter = float(full_filter)
    relaxed_filter = float(relaxed_filter)

    # Basic cleanup
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)

    # QC annotations
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB(?!P)", regex=True)

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True)

    qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

    sc.settings.figdir = violin_dir
    sc.pl.violin(adata, qc_cols, save=f"_{subproject_name}_QC.png", jitter=0.4, multi_panel=True)

    sc.settings.figdir = scatter_dir
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save=f"_{subproject_name}_QC.png")

    # Adaptive filtering (important for tumor heterogeneity)
    if adata.n_obs > 50:
        cell_mask = (
            (adata.obs["n_genes_by_counts"] >= adata.obs["n_genes_by_counts"].quantile(full_filter)) &
            (adata.obs["total_counts"] >= adata.obs["total_counts"].quantile(full_filter)) &
            (adata.obs["pct_counts_mt"] <= adata.obs["pct_counts_mt"].quantile(1 - full_filter))
        )
    else:
        cell_mask = (
            (adata.obs["n_genes_by_counts"] >= adata.obs["n_genes_by_counts"].quantile(relaxed_filter)) &
            (adata.obs["total_counts"] >= adata.obs["total_counts"].quantile(relaxed_filter)) &
            (adata.obs["pct_counts_mt"] <= adata.obs["pct_counts_mt"].quantile(1 - relaxed_filter))
        )

    adata = adata[cell_mask].copy()

    # Preserve raw counts
    if "raw_counts" in adata.layers:
        adata.layers["counts"] = adata.layers["raw_counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    # Gene filtering (keep important genes)
    genes_expressed_mask = sc.pp.filter_genes(adata, min_cells=3, inplace=False)[0]
    important_set = {g.upper() for g in important_genes}
    important_genes_mask = adata.var_names.str.upper().isin(important_set)

    final_gene_mask = genes_expressed_mask | important_genes_mask
    adata = adata[:, final_gene_mask].copy()

    # Normalize AFTER filtering (critical difference from embryo pipeline)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store clean reference
    adata.raw = adata.copy()

    adata.write(adata_filt_path)

    return adata

def concat_cancer_adata(concat_adatas, concat_keys, concat_path): 
    """After filtering (QC, genes, cells, keep important genes) but before 
    Normalize/Log; HVG; PCA;Neighbors;UMAP;clustering"""
    
    
    adata_concat = ad.concat(concat_adatas, axis=0, join="outer", label="dataset", keys=concat_keys, fill_value=0)
    adata_concat.write(concat_path)
    return adata_concat

def run_qc_ovarian(project_name, adata_filtered, adata_qc_path, important_genes, qc_save_dir, cell_cycle_genes_file_path):

    pca_dir = os.path.join(qc_save_dir, "pca")
    hvg_dir = os.path.join(qc_save_dir, "hvg")
    os.makedirs(pca_dir, exist_ok=True)
    os.makedirs(hvg_dir, exist_ok=True)

    n_before = adata_filtered.n_obs
    adata = adata_filtered.copy()

    # Ensure batch key exists (for concat)
    if "sample_id" not in adata.obs:
        adata.obs["sample_id"] = adata.obs["dataset"].astype(str)

    # Remove doublets FIRST
    if "predicted_doublet" in adata.obs:
        adata = adata[~adata.obs["predicted_doublet"]].copy()

    n_after = adata.n_obs
    print(f"Removed {n_before - n_after} doublets in project {project_name}")

    # Restore raw counts into X + ensure counts layer exists
    if "raw_counts" in adata.layers:
        adata.X = adata.layers["raw_counts"].copy()
        adata.layers["counts"] = adata.layers["raw_counts"].copy()
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    # Remove dead genes (important after concat)
    sc.pp.filter_genes(adata, min_cells=3)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    sc.settings.figdir = hvg_dir
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample_id", flavor="seurat")
    sc.pl.highly_variable_genes(adata, save=f"_{project_name}_HVGs.png")

    # Subset to HVGs (IMPORTANT)
    adata = adata[:, adata.var["highly_variable"]].copy()

    adata.raw = adata.copy()

    # Cell cycle scoring
    cell_cycle_scoring(project_name, adata, cell_cycle_genes_file_path, qc_save_dir)

    sc.pp.scale(adata, max_value=10)

    max_pcs = max(2, min(adata.n_obs, adata.n_vars) - 1)
    n_pcs = min(30, max_pcs)
    sc.tl.pca(adata, n_comps=n_pcs)

    sc.settings.figdir = pca_dir
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save=f"_{project_name}_PCA_variance_ratio.png")

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)

    for res in [0.5, 1.0, 2.0, 3.0]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:3.1f}", resolution=res, flavor="igraph")

    sc.tl.umap(adata)

    adata.write(adata_qc_path)

    return adata