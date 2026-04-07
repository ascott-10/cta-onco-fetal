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
from codes.import_raw_data import import_raw_data_fetal_gonad, import_raw_data_embryos_mixed, import_raw_data_10x_subprojects
from codes.filter_qc_data import filter_data,filter_data_embryos_mixed, run_qc, run_qc_embryos_mixed
from codes.cell_annotation import cell_annotation_wrapper, annotate_obs_columns, adata_split_annotate
from codes.visual_gene_expression import make_umaps, cta_genes_expression, cta_genes_expression_all_samples
from codes.differential_gene_expression import get_ranked_genes, markers_leiden_embryos_mixed, rank_genes_for_publication

def pre_process_project_setup(project_name):
    cfg = set_up_project_config(project_name)
    paths, project_cfg, global_cfg = cfg["paths"], cfg["project_cfg"], cfg["global_cfg"]
    

    # Set up samples
    
    samples_to_process = project_cfg["sample_list"]
    

    print("# Samples to Process: ", len(samples_to_process))

    return paths, project_cfg, global_cfg, samples_to_process

def run_process_pipeline(project_name, paths, project_cfg, global_cfg, samples_to_process):

    important_genes = global_cfg["important_genes"]
    cell_cycle_genes_file_path = CELL_CYCLE_GENES_FILE_PATH
    cols_to_plot = ["sample_id", "sex", "stage_label"]
    cols_to_plot_titles = ["Sample", "Sex", "Age (Gestation Weeks)"]
    final_leiden_res = global_cfg["final_leiden_res"]
    json_annotations_path = global_cfg["json_annotations_path"]

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
                adata_filtered = filter_data_embryos_mixed(project_name, adata_init, adata_filt_path, important_genes, cell_cycle_genes_file_path, paths["QC_SAVE_DIR"])   
            adata_qc = run_qc_embryos_mixed(project_name, adata_filtered, adata_qc_path, important_genes, paths["QC_SAVE_DIR"])
        
            
            
        adata_female, adata_male = adata_split_annotate(project_name, final_leiden_res, adata_female_path, adata_male_path, adata_qc_path, json_annotations_path)
        for sex, adata_current, adata_path in [
            ("female", adata_female, adata_female_path),
            ("male", adata_male, adata_male_path),]:
            subproject_name = f"{project_name}_{sex}"
            adata_annotated = annotate_obs_columns(subproject_name, adata_current, cols_to_plot, cols_to_plot_titles, adata_path, paths["FIGURES_DIR"], palette = "tab10")
            top_markers_path = os.path.join(paths["TABLES_DIR"], "deg_analysis", f"leiden_{final_leiden_res}_{sex}_markers.csv")
            markers_leiden_embryos_mixed(adata_annotated, subproject_name, top_markers_path, final_leiden_res, paths["FIGURES_DIR"])