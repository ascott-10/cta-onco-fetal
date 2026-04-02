# python_scripts_common/main.py

from config import *


import os

import scanpy as sc
import anndata as ad

os.chdir("/data/scottaa/cta_onco_fetal")

# main.py
def main():
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

if __name__ == "__main__":
    main()



# source myconda; conda activate rnaseq-pipe
# python python_common_scripts/concat_adata.py