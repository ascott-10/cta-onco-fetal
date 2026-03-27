# python_scripts_common/codes/import_raw_data.py

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

#### Functions ####

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

def import_raw_data_csv(project_name, sample_id, original_data_dir,adata_init_path, sample_meta_df):
    """    ["embryos_mixed"] """
    
    if project_name in ["embryos_mixed"]:
        raw_counts_path = os.path.join(original_data_dir, f"{sample_id}_concat_gene_expression.txt")
    df = pd.read_csv(raw_counts_path, index_col = 0).T
    
    adata = sc.AnnData(df)
    adata.var_names_make_unique()
 
    adata.layers["raw_counts"] = adata.X.copy()
    adata.raw = adata.copy()

    adata.obs["sample_id"] = sample_id
    meta_row = sample_meta_df.set_index("sample_id").loc[sample_id]

    for col in sample_meta_df.columns:
        if col == "sample_id":
            continue
        adata.obs[col] = meta_row[col]
    
    adata.write(adata_init_path)

    return adata
    
  

def import_raw_data_10x(subproject, original_data_dir,adata_init_path):

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