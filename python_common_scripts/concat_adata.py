# python_scripts_common/main.py

from config import *


import os

import scanpy as sc
import anndata as ad

os.chdir("/data/scottaa/cta_onco_fetal")

# main.py


def main_embryos():
    adata_embryos = sc.read_h5ad("datasets/embryos_mixed/working_adata/embryos_mixed_concat_annotated_cellmarker.h5ad")
    adata_fetal = sc.read_h5ad("datasets/fetal_gonad/working_adata/fetal_gonad_concat_annotated_cellmarker.h5ad")

    adatas = [adata_embryos, adata_fetal]
    adata_key_list= ["embryos_mixed", "fetal_gonad"]


    adata_concat = ad.concat(adatas,
    axis=0,          # 0 = concatenate cells, 1 = concatenate genes
    join="outer",    # 'outer' keeps all genes, 'inner' keeps only shared genes
    label="dataset",   # Adds a column in .obs to indicate source
    keys=adata_key_list,  # Names for each batch
    fill_value=0     # Fill missing values with 0
    )
    adata_concat.write("datasets/all_fetal_concat_annotated_cellmarker.h5ad")
    print(adata_concat)

def main_cancer():
    working_adata_dir = "datasets/ovarian_cancer_ccca/working_adata"
    raw_data_dir = "datasets/ovarian_cancer_ccca/raw_data"
    concat_adatas = []
    concat_keys = []
    concat_path = os.path.join(working_adata_dir, "all_ovarian_cancer_concat_before_qc.h5ad")
    for adata_path in os.listdir(raw_data_dir):
        if adata_path.endswith("after_hvg.h5ad"):
            adata_name = adata_path.split("_")[0]
            full_adata_path = os.path.join(raw_data_dir, adata_path)
            adata = sc.read_h5ad(full_adata_path)
            concat_adatas.append(adata)
            concat_keys.append(adata_name)

    adata_concat = ad.concat(concat_adatas,axis=0, join="outer", label="dataset", keys=concat_keys, fill_value=0)
    adata_concat.write(concat_path)
    print(adata_concat)

if __name__ == "__main__":
    #main_embryos()
    main_cancer()



# source myconda; conda activate rnaseq-pipe
# python python_common_scripts/concat_adata.py