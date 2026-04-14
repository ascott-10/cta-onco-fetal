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

    if os.path.exists(df_all_save_path):
        df_all = pd.read_csv(df_all_save_path)

    else:
        for file in os.listdir(original_original_data_dir):
            if file.endswith(".txt"):
                file_name = os.path.join(original_original_data_dir, file)
                df = pd.read_csv(file_name, sep='\t')
                df = df.set_index("Gene")
                original_file_df_list.append(df)

        df_all = pd.concat(original_file_df_list, axis=1).fillna(0)
        df_all.reset_index().to_csv(df_all_save_path, index=False)

    df_all = df_all.set_index("Gene")

    adata = sc.AnnData(df_all.T)
    adata.var_names_make_unique()

    adata.layers["raw_counts"] = adata.X.copy()
    adata.raw = adata.copy()

    adata.obs["sample_id"] = adata.obs.index.str.split("_embryo").str[0]
    sample_meta = pd.read_csv(sample_meta_path)
    for col in adata.obs.columns:
        if adata.obs[col].dtype == "object":
            try:
                adata.obs[col] = adata.obs[col].astype("string")
            except:
                print(f"Could not convert {col}")
    sample_meta = sample_meta.set_index("sample_id")

    adata.obs = adata.obs.join(sample_meta, on ="sample_id")

    adata.write(adata_init_path)

    return adata

def import_raw_data_10x_ovarian_cancer(subproject, original_data_dir,adata_init_path, sample_meta_path):

    """["ovarian_cancer_ccca"]
    ['Nath2021','Olalekan2021','Olbrecht2021','Qian2020','Regner2021','Shih2018','Tang-Huau2018','Zhang2019', 'Zhang2022']"""

    
    matrix_path = os.path.join(original_data_dir, "Exp_data_UMIcounts.mtx")
    barcodes_path = os.path.join(original_data_dir,  "Cells.csv")
    genes_path = os.path.join(original_data_dir,  "Genes.txt")

    from scipy.io import mmread

    with open(matrix_path, "rb") as f:
        # skip until we hit the real header
        for line in f:
            if line.startswith(b"%%MatrixMarket"):
                break
        # now read from correct position
    X = mmread(f).tocsr().T
   

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
    cols_to_add = [c for c in sample_meta_df.columns if c not in adata.obs.columns or c == "sample_id"]

    adata.obs["sample_id"] = adata.obs["sample"].astype(str)
    sample_meta_df["sample_id"] = sample_meta_df["sample_id"].astype(str)


    adata.obs = adata.obs.reset_index().merge(sample_meta_df[cols_to_add], on="sample_id", how="left").set_index("cell_name")
    
    
    adata.write(adata_init_path)

    return adata




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

  

