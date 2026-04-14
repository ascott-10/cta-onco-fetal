# python_scripts_common/main_process_data.py

import os
import glob
import sys
import pandas as pd
import scanpy as sc
import json
import scanpy as sc

# Import from your config and helper codes
from config import *
from codes.import_raw_data import import_raw_data_embryos_mixed, import_raw_data_10x_ovarian_cancer
from codes.filter_qc_data import filter_data_embryos_mixed, filter_data_ovarian, run_qc_embryos_mixed, run_qc_ovarian, concat_cancer_adata
from codes.cell_annotation import adata_split, run_full_annotation_pipeline
from codes.visual_gene_expression import show_cta_genes
#from codes.visual_gene_expression import make_umaps, cta_genes_expression, cta_genes_expression_all_samples
#from codes.differential_gene_expression import get_ranked_genes, markers_leiden_embryos_mixed, rank_genes_for_publication

def pre_process_project_setup(project_name):
    cfg = set_up_project_config(project_name)
    paths, project_cfg, global_cfg = cfg["paths"], cfg["project_cfg"], cfg["global_cfg"]
    

    # Set up samples
    
    samples_to_process = project_cfg["sample_list"]
    

    print("# Samples to Process: ", len(samples_to_process))

    return paths, project_cfg, global_cfg, samples_to_process

def run_process_pipeline(project_name, paths, project_cfg, global_cfg, samples_to_process):

    important_genes = global_cfg["important_genes"]
    leiden_res_list = global_cfg["leiden_res_list"]
    json_annotations_path = global_cfg["json_annotations_path"]
    palette_path = global_cfg["json_palette_path"]

    cols_to_plot = project_cfg["plot"]["obs_cols"]
    cols_to_plot_titles = project_cfg["plot"]["obs_titles"]

    cell_cycle_genes_file_path = global_cfg["cell_cycle_genes_path"]

    if project_name in ["embryos_mixed"]:
        original_original_data_dir = os.path.join(paths["ORIGINAL_ORIGINAL_DATA_DIR"])
        original_data_dir = os.path.join(paths["ORIGINAL_DATA_DIR"])
        sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])

        adata_init_path = os.path.join(paths["RAW_DATA_DIR"], f"{project_name}_init.h5ad")
        adata_filt_path = os.path.join(paths["RAW_DATA_DIR"], f"{project_name}_before_qc.h5ad")
        adata_qc_path = os.path.join(paths["RAW_DATA_DIR"], f"{project_name}_mixed_after_hvg.h5ad")
        #adata_annotated_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_mixed_processed.h5ad")
        adata_male_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_male_processed.h5ad")
        adata_female_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_female_processed.h5ad")

        adata_male_annotated_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_male_processed_annotated.h5ad")
        adata_female_annotated_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_female_processed_annotated.h5ad")

        
        if os.path.exists(adata_qc_path):
            adata_qc = sc.read_h5ad(adata_qc_path)
        else:
            if os.path.exists(adata_filt_path):
                adata_filtered = sc.read_h5ad(adata_filt_path)
            else:
                if os.path.exists(adata_init_path):
                    adata_init = sc.read_h5ad(adata_init_path)
                else:
                    adata_init = import_raw_data_embryos_mixed(project_name, original_original_data_dir, original_data_dir, adata_init_path, sample_meta_path)
                adata_filtered = filter_data_embryos_mixed(project_name, adata_init, adata_filt_path, important_genes, paths["QC_SAVE_DIR"])
            adata_qc = run_qc_embryos_mixed(project_name, adata_filtered, adata_qc_path, important_genes, paths["QC_SAVE_DIR"],cell_cycle_genes_file_path)
         
        
        
        adata_female, adata_male = adata_split(adata_female_path, adata_male_path, adata_female_annotated_path,adata_male_annotated_path, adata_qc_path)
        
        for sex, adata_current,adata_annotated_path in [("female", adata_female, adata_female_annotated_path),
                                                ("male", adata_male, adata_male_annotated_path),]:
            subproject_name = f"{project_name}_{sex}"
            
            if os.path.exists(adata_annotated_path):
                adata_annotated = sc.read_h5ad(adata_annotated_path)      
                print(adata_annotated.obs.columns)
            else:
                
                adata_annotated = run_full_annotation_pipeline(subproject_name,adata_current.copy(), adata_annotated_path, leiden_res_list,json_annotations_path, palette_path,cols_to_plot,cols_to_plot_titles,tables_dir=paths["TABLES_DIR"],figures_dir=paths["FIGURES_DIR"],palette="tab10")
                print(adata_annotated.obs.columns)
            show_cta_genes(subproject_name, adata_annotated, figures_dir=paths["FIGURES_DIR"], tables_dir=paths["TABLES_DIR"],json_path = global_cfg["json_annotations_path"], groupby = "celltype_leiden_res_1.0")
            
    elif project_name in ["ovarian_cancer_ccca"]:
        samples_to_process = ['Nath2021','Olalekan2021','Olbrecht2021','Qian2020','Regner2021','Shih2018','Tang-Huau2018','Zhang2019', 'Zhang2022']
        leiden_res_list=[0.5]
        concat_filt_path = os.path.join(paths["WORKING_ADATA_DIR"], "all_ovarian_cancer_concat_before_qc.h5ad")
        concat_qc_path = os.path.join(paths["WORKING_ADATA_DIR"], "all_ovarian_cancer_concat_after_hvg.h5ad")
        concat_annotated_path = os.path.join(paths["WORKING_ADATA_DIR"], "all_ovarian_cancer_concat_annotated.h5ad")
        
        concat_adatas = []
        concat_keys = []

        # Check if the concat adata was already made; if so, load the most annotated version
        if os.path.exists(concat_qc_path):
            adata_concat_qc = sc.read_h5ad(concat_qc_path)
        else:
            # If only a raw version was made, load it
            if os.path.exists(concat_filt_path):
                adata_concat_filt = sc.read_h5ad(concat_filt_path)
                # If not, then must build the concat
            else:

                for subproject in samples_to_process:

                    original_data_dir = os.path.join(paths["ORIGINAL_DATA_DIR"], f"Data_{subproject}_Ovarian")
                    sample_meta_path = os.path.join(original_data_dir, "Samples copy.csv")

                    adata_init_path = os.path.join(paths["RAW_DATA_DIR"], f"{subproject}_init.h5ad")
                    adata_filt_path = os.path.join(paths["RAW_DATA_DIR"], f"{subproject}_before_qc.h5ad")

                    # Check if individual first pass qc exists; load it
                    if os.path.exists(adata_filt_path):
                        print(f"reading {subproject}")
                        adata_filtered = sc.read_h5ad(adata_filt_path)
                    else:
                        # Otherwise, check if initial processed exists; load it
                        if os.path.exists(adata_init_path):
                            adata_init = sc.read_h5ad(adata_init_path)
                        else:
                            # Else import raw data
                            adata_init = import_raw_data_10x_ovarian_cancer(subproject, original_data_dir, adata_init_path, sample_meta_path)
                        # Then filter 
                        adata_filtered = filter_data_ovarian(subproject, adata_init, adata_filt_path, important_genes, paths["QC_SAVE_DIR"], full_filter=0.05, relaxed_filter=0.02)

                        
                    # Append either the read-in filtered, or the created filtered
                    concat_adatas.append(adata_filtered)
                    concat_keys.append(subproject)
                # Finally, use the lists to create concat adatas
                print("concat keys list: ",concat_keys)
                adata_concat_filt = concat_cancer_adata(concat_adatas, concat_keys, concat_filt_path)
                # Then do qc on the final concat
                adata_concat_qc = run_qc_ovarian(project_name, adata_concat_filt, concat_qc_path, important_genes, paths["QC_SAVE_DIR"], cell_cycle_genes_file_path)
        
        # Finally run annotation pipeline on 
        adata_concat_annotated = run_full_annotation_pipeline(project_name,adata_concat_qc.copy(), concat_annotated_path, leiden_res_list,json_annotations_path, palette_path,cols_to_plot,cols_to_plot_titles,tables_dir=paths["TABLES_DIR"],figures_dir=paths["FIGURES_DIR"],palette="tab10")

       

        fetal_genes = pd.read_csv("datasets/embryos_mixed/results/tables/cta_analysis/embryos_mixed_female_cta_genes.csv")["gene"].tolist()
        show_cta_genes(project_name, adata_concat_annotated, figures_dir=paths["FIGURES_DIR"], tables_dir=paths["TABLES_DIR"],json_path = global_cfg["json_annotations_path"], groupby = "celltype_leiden_res_0.5", cta_genes_fixed=fetal_genes)
    