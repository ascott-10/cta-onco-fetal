# python_scripts_common/codes/import_raw_data.py

#### Libraries ###
import os
import glob
from io import StringIO
import pandas as pd
import numpy as np
import scanpy as sc

import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread


import matplotlib.pyplot as plt

#### Custom ####


from config import *

#### Functions ####

def import_raw_data_embryos_mixed(project_name, original_original_data_dir, original_data_dir, adata_init_path, sample_meta_path):
    # Concat all text files
    original_file_df_list = []
    df_all_save_path = os.path.join(original_data_dir, f"{project_name}_gene_expression.csv")

    for file in os.listdir(original_original_data_dir):
        if file.endswith(".txt"):
            file_name = os.path.join(original_original_data_dir, file)
            df = pd.read_csv(file_name, sep='\t')
            df = df.set_index("Gene")
            original_file_df_list.append(df)

    df_all = pd.concat(original_file_df_list, axis=1).fillna(0)
    df_all.reset_index().to_csv(df_all_save_path, index=False)

    df = df.set_index("Gene")

    adata = sc.AnnData(df.T)
    adata.var_names_make_unique()

    adata.layers["raw_counts"] = adata.X.copy()
    adata.raw = adata.copy()

    adata.obs["sample_id"] = adata.obs.index.str.split("_embryo").str[0]
    sample_meta = pd.read_csv(sample_meta_path)
    sample_meta = sample_meta.set_index("sample_id")

    adata.obs = adata.obs.join(sample_meta, on ="sample_id")

    adata.write(adata_init_path)


def import_raw_data_fetal_gonad(sample_id, gse_id, original_data_dir,adata_init_path, sample_meta_df):

    """["fetal_gonad]"""
    
    if gse_id in ['GSM5506062_G1', 'GSM5506063_G2','GSM5506064_G3', 'GSM5506065_G4', 'GSM5704349_mesonephros']:
        
        adata = sc.read_10x_h5(os.path.join(original_data_dir, f"{gse_id}.h5"))
    
    elif gse_id in ['GSM6703999_mesonephros_F','GSM6704000_G5_A','GSM6704001_G5_B','GSM6704002_G6_A','GSM6704003_G6_B']:
    
        matrix_path = os.path.join(original_data_dir, gse_id, "matrix.mtx")
        genes_path = os.path.join(original_data_dir, gse_id, "features.tsv")
        barcodes_path = os.path.join(original_data_dir, gse_id, "barcodes.tsv")
        
        X = mmread(matrix_path).tocsr().T

        genes = pd.read_csv(genes_path, header=None, sep="\t")
        barcodes = pd.read_csv(barcodes_path, header=None)

        adata = sc.AnnData(X)

        adata.var_names = genes[1].astype(str).values
        adata.obs_names = barcodes[0].astype(str).values

    print("var_names",adata.var_names[:10])
    print("obs_names", adata.obs_names[:10])

    adata.var_names_make_unique()

    adata.layers["raw_counts"] = adata.X.copy()
    adata.raw = adata.copy()

    adata.obs["sample_id"] = sample_id
    print(adata.obs)
    sample_meta_df.set_index("sample_id", drop = False, inplace=True)
    print(sample_meta_df)
    row = sample_meta_df.loc[sample_id]
    for col, value in row.items():
        adata.obs[col] = value
    
    adata.write(adata_init_path)

    return adata

  

def import_raw_data_10x_subprojects(subproject, original_data_dir,adata_init_path):

    """["ovarian_cancer_ccca"]"""

    if subproject in ["Izar2020"]:
        matrix_path = os.path.join(original_data_dir, f"Data_{subproject}_Ovarian", "Exp_data_TPM.mtx")
    else:
        matrix_path = os.path.join(original_data_dir, f"Data_{subproject}_Ovarian", "Exp_data_UMIcounts.mtx")
    
    barcodes_path = os.path.join(original_data_dir, f"Data_{subproject}_Ovarian", "Cells.csv")
    genes_path = os.path.join(original_data_dir, f"Data_{subproject}_Ovarian", "Genes.txt")
    sample_meta_path = os.path.join(original_data_dir, f"Data_{subproject}_Ovarian", "Samples.csv")

    X = mmread(matrix_path).tocsr().T

    genes = pd.read_csv(genes_path, header=None)[0].astype(str).str.replace('"', '', regex=False)
    barcodes = pd.read_csv(barcodes_path, sep=None, engine="python")

    adata = sc.AnnData(X)
    adata.var_names = genes.values
    adata.obs_names = barcodes["cell_name"].astype(str).values
    adata.obs = barcodes.set_index("cell_name")

    adata.layers["raw_counts"] = adata.X.copy()
    adata.raw = adata.copy()

    sample_meta_df = pd.read_csv(sample_meta_path)
    sample_meta_df = sample_meta_df.dropna(axis=1, how="all")
    sample_meta_df = sample_meta_df.astype(str)

    
    # columns to add (exclude overlaps except 'sample')
    cols_to_add = [c for c in sample_meta_df.columns if c not in adata.obs.columns or c == "sample"]

    adata.obs["sample"] = adata.obs["sample"].astype(str)
    sample_meta_df["sample"] = sample_meta_df["sample"].astype(str)


    adata.obs = adata.obs.reset_index().merge(sample_meta_df[cols_to_add], on="sample", how="left").set_index("cell_name")
    adata.obs = adata.obs.rename(columns={"sample": "sample_id"})
    
    adata.write(adata_init_path)

    return adata


def import_raw_data_10x(project_name, sample_id, gse_id, original_data_dir,adata_init_path, sample_meta_df):

    """["subtype_evolution", "hgsoc_subtype_define", "hgsoc_tissue_architecture", "gyne_malignant"]"""

    
    if project_name in ["subtype_evolution"]:
        
        matrix_path = os.path.join(original_data_dir, f"{gse_id}_{sample_id}_matrix.mtx")
        genes_path = os.path.join(original_data_dir, f"{gse_id}_{sample_id}_genes.tsv")
        barcodes_path = os.path.join(original_data_dir, f"{gse_id}_{sample_id}_barcodes.tsv")
        
    elif project_name in ["hgsoc_subtype_define"]:

        matrix_path = os.path.join(original_data_dir, f"{gse_id}_matrix_{sample_id}.mtx")
        genes_path = os.path.join(original_data_dir, f"{gse_id}_features_{sample_id}.tsv")
        barcodes_path = os.path.join(original_data_dir, f"{gse_id}_barcodes_{sample_id}.tsv")
        
    elif project_name in ["hgsoc_tissue_architecture"]:

        matrix_path = os.path.join(original_data_dir, f"{gse_id}_{sample_id}.matrix.mtx")
        genes_path = os.path.join(original_data_dir, f"{gse_id}_{sample_id}.genes.tsv")
        barcodes_path = os.path.join(original_data_dir, f"{gse_id}_{sample_id}.barcodes.tsv")

    elif project_name in ["gyne_malignant"]:

        matrix_path = os.path.join(original_data_dir, f"{gse_id}_matrix-{sample_id}.mtx")
        genes_path = os.path.join(original_data_dir, f"{gse_id}_features-{sample_id}.tsv")
        barcodes_path = os.path.join(original_data_dir, f"{gse_id}_barcodes-{sample_id}.tsv")
        
        
    genes = pd.read_csv(genes_path, sep='\t', header=None)
    barcodes = pd.read_csv(barcodes_path, sep='\t', header=None)
    barcodes.columns = ["cell_id"]

    X = mmread(matrix_path).tocsr().T

    # shape checks
    assert X.shape[1] == genes.shape[0], f"Mismatch: {X.shape[1]} genes in matrix vs {genes.shape[0]} in genes.tsv"
    assert X.shape[0] == barcodes.shape[0], f"Mismatch: {X.shape[0]} cells in matrix vs {barcodes.shape[0]} in barcodes.tsv"

    adata = sc.AnnData(X)

    
    adata.var["ensembl_id"] = genes[0].astype(str).values
    adata.var["gene_symbol"] = genes[1].astype(str).values
    adata.var["gene_type"] = genes[2].astype(str).values

    adata.var_names = adata.var["gene_symbol"]
    adata.var_names.name = None

    adata.var.drop(columns="gene_symbol", inplace=True)
    
    adata.var_names_make_unique()

    adata.obs_names = barcodes["cell_id"].astype(str).values
    adata.obs = barcodes.set_index("cell_id")

    adata.layers["raw_counts"] = adata.X.copy()
    adata.raw = adata.copy()

    adata.obs["sample_id"] = sample_id
    
    sample_meta_df = sample_meta_df.copy()
    sample_meta_df.set_index("sample_id", drop=False, inplace=True)

    row = sample_meta_df.loc[sample_id]
    for col, value in row.items():
        adata.obs[col] = value
    
    adata.write(adata_init_path)

    return adata