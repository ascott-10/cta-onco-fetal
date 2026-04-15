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
from codes.import_raw_data import import_raw_data_embryos_mixed, import_raw_data_10x_ensembl, import_raw_data_h5_gsm, concat_adata
from codes.filter_qc_data import filter_data_embryos_mixed, run_qc_embryos_mixed, filter_data, run_qc
from codes.cell_annotation import adata_split, run_full_annotation_pipeline
from codes.visual_gene_expression import show_cta_genes, build_palettes, plot_umaps
from codes.differential_gene_expression import get_ranked_genes

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
    gene_map_file_path = global_cfg["gene_map_file_path"]

    cols_to_plot = project_cfg["plot"]["obs_cols"]
    cols_to_plot_titles = project_cfg["plot"]["obs_titles"]

    cell_cycle_genes_file_path = global_cfg["cell_cycle_genes_path"]
    cta_genes_file_path = global_cfg["cta_genes_path"]

    if os.path.exists(palette_path):
        try:
            with open(palette_path, "r") as f:
                existing_palettes = json.load(f)
        except Exception:
            existing_palettes = {}
    else:
        existing_palettes = {}

    if project_name in ["fetal_gonad"]:

        original_original_data_dir = os.path.join(paths["ORIGINAL_ORIGINAL_DATA_DIR"])
        original_data_dir = os.path.join(paths["ORIGINAL_DATA_DIR"])
        sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])

        gse_list = project_cfg["gse_list"]
        sample_list = project_cfg["sample_list"]
        group_list = project_cfg["group_list"]

        adata_male = []
        adata_male_key = []
        adata_female = []
        adata_female_key = []
        adata_mixed = []
        adata_mixed_key = []

        for gse_id, sample_id in zip(gse_list, sample_list):

            adata_init_path = os.path.join(original_data_dir, f"{sample_id}_init.h5ad")

            if os.path.exists(adata_init_path):
                adata_init = sc.read_h5ad(adata_init_path)
            else:
                adata_init = import_raw_data_h5_gsm(project_name, gse_id, original_original_data_dir, adata_init_path, sample_meta_path)

            if sample_id in ['G1_1_Ovary', 'G2_2_Ovary','G_2_F_Mesonephros','G5_A_2_Ovary','G5_B_2_Ovary']:
                adata_female.append(adata_init)
                adata_female_key.append(sample_id)
     
            elif sample_id in ['G3_1_Testis','G4_2_Testis','G_1_2_M_Mesonephros']:
                adata_male.append(adata_init)
                adata_male_key.append(sample_id)

            elif sample_id in ['G6_A_Mixed','G6_B_Mixed']:
                adata_mixed.append(adata_init)
                adata_mixed_key.append(sample_id)


        for group, adata_list, adata_key_list in [("Female", adata_female, adata_female_key),
                                                ("Male", adata_male, adata_male_key),
                                                ("Mixed", adata_mixed, adata_mixed_key),]:
            
            subproject_name = f"{project_name}_{group}"
            
            adata_concat_init_path = os.path.join(paths["RAW_DATA_DIR"], f"{subproject_name}_init.h5ad") 
            adata_filt_path = os.path.join(paths["RAW_DATA_DIR"], f"{subproject_name}_before_qc.h5ad")
            adata_qc_path = os.path.join(paths["RAW_DATA_DIR"], f"{subproject_name}_after_hvg.h5ad")
            adata_ranked_path = os.path.join(paths["RAW_DATA_DIR"], f"{subproject_name}_ranked.h5ad")
            adata_annotated_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{subproject_name}_annotated.h5ad")


            if os.path.exists(adata_qc_path):
                adata_qc = sc.read_h5ad(adata_qc_path)
            else:
                if os.path.exists(adata_filt_path):
                    adata_filtered = sc.read_h5ad(adata_filt_path)         
                else:
                    if os.path.exists(adata_concat_init_path):
                        adata_concat_init = sc.read_h5ad(adata_concat_init_path)
                    else:
                        adata_concat_init = concat_adata(adata_list, adata_key_list, adata_concat_init_path)
                    adata_filtered = filter_data(subproject_name, adata_concat_init, adata_filt_path, important_genes, paths["QC_SAVE_DIR"], full_filter=0.05, relaxed_filter=0.02)
            
                adata_qc = run_qc(subproject_name, adata_filtered, adata_qc_path, paths["QC_SAVE_DIR"], cell_cycle_genes_file_path)

        
            adata_ranked, leiden_cols, leiden_cols_titles, celltype_cols = get_ranked_genes(adata_qc, adata_ranked_path, subproject_name, json_annotations_path, leiden_res_list, paths["TABLES_DIR"],figures_dir=paths["FIGURES_DIR"], always_rank=True)

            all_cols = cols_to_plot + leiden_cols + celltype_cols
            plot_cols = cols_to_plot + celltype_cols
            plot_titles = (cols_to_plot_titles + [f"Predicted Cell Type \n Leiden Resolution {c.replace('celltype_leiden_res_', '')}" for c in celltype_cols])
        
            adata_annotated, existing_palettes = build_palettes(adata_ranked, adata_annotated_path, all_cols,existing_palettes)

    
            plot_umaps(adata_annotated, adata_annotated_path, subproject_name,plot_cols, plot_titles, leiden_cols,leiden_cols_titles,figures_dir=paths["FIGURES_DIR"])

           

            cta_genes_list = show_cta_genes(subproject_name, adata_annotated, figures_dir=paths["FIGURES_DIR"], tables_dir=paths["TABLES_DIR"], cta_genes_file_path = cta_genes_file_path, groupby = "celltype_leiden_res_1.0", cta_genes_fixed=None)
    
            adata_annotated.write(adata_annotated_path)
        
        with open(palette_path, "w") as f:
            json.dump(existing_palettes, f, indent=2)
        
            


        




    if project_name in ["mtab_tumors"]:

        original_data_dir = os.path.join(paths["ORIGINAL_DATA_DIR"])
        sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])

        adata_init_path = os.path.join(paths["RAW_DATA_DIR"], f"{project_name}_init.h5ad")
        adata_filt_path = os.path.join(paths["RAW_DATA_DIR"], f"{project_name}_before_qc.h5ad")
        adata_qc_path = os.path.join(paths["RAW_DATA_DIR"], f"{project_name}_after_hvg.h5ad")
        adata_ranked_path = os.path.join(paths["RAW_DATA_DIR"], f"{project_name}_ranked.h5ad")
        adata_annotated_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_annotated.h5ad")

        
        if os.path.exists(adata_qc_path):
            adata_qc = sc.read_h5ad(adata_qc_path)
        else:
            if os.path.exists(adata_filt_path):
                adata_filtered = sc.read_h5ad(adata_filt_path)
            else:
                if os.path.exists(adata_init_path):
                    adata_init = sc.read_h5ad(adata_init_path)
                else:
                    adata_init = import_raw_data_10x_ensembl(project_name, original_data_dir, adata_init_path, sample_meta_path, gene_map_file_path)
                adata_filtered = filter_data(project_name, adata_init, adata_filt_path, important_genes, paths["QC_SAVE_DIR"], full_filter=0.05, relaxed_filter=0.02)
            adata_qc = run_qc(project_name, adata_filtered, adata_qc_path, paths["QC_SAVE_DIR"], cell_cycle_genes_file_path)

        

        
        adata_ranked, leiden_cols, leiden_cols_titles, celltype_cols = get_ranked_genes(adata_qc, adata_ranked_path, project_name, json_annotations_path, leiden_res_list, paths["TABLES_DIR"],figures_dir=paths["FIGURES_DIR"], always_rank=True)

        all_cols = cols_to_plot + leiden_cols + celltype_cols
        plot_cols = cols_to_plot + celltype_cols
        plot_titles = (cols_to_plot_titles + [f"Predicted Cell Type \n Leiden Resolution {c.replace('celltype_leiden_res_', '')}" for c in celltype_cols])
        
        adata_annotated, all_palettes = build_palettes(adata_ranked, adata_annotated_path, all_cols,existing_palettes)

    
        plot_umaps(adata_annotated, adata_annotated_path, project_name,plot_cols, plot_titles, leiden_cols,leiden_cols_titles,figures_dir=paths["FIGURES_DIR"])

        with open(palette_path, "w") as f:
            json.dump(all_palettes, f, indent=2)
        
        adata_annotated.write(adata_annotated_path)

        cta_genes_list = show_cta_genes(project_name, adata_annotated, figures_dir=paths["FIGURES_DIR"], tables_dir=paths["TABLES_DIR"], cta_genes_file_path = cta_genes_file_path, groupby = "celltype_leiden_res_1.0", cta_genes_fixed=None)
    
        adata_annotated.write(adata_annotated_path)
        

         
        

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
            
    
      