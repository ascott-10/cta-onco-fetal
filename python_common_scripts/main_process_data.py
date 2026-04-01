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
from codes.import_raw_data import import_raw_data_fetal_gonad, import_raw_data_10x, import_raw_data_10x_subprojects, import_raw_data_csv
from codes.filter_qc_data import filter_data, run_qc
from codes.cell_annotation import cell_annotation_wrapper
from codes.visual_gene_expression import make_umaps, cta_genes_expression, cta_genes_expression_all_samples
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
    
        
    
    
    for sample_id in samples_to_process:
        print("Processing sample", sample_id)

        # Set sample meta path
        if project_name in ["subtype_evolution", "hgsoc_subtype_define", "hgsoc_tissue_architecture", "gyne_malignant"]:
            # if '\t' separated
            sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])
            sample_meta_df = pd.read_csv(sample_meta_path, usecols=project_cfg["sample_meta_cols"], sep = '\t')
        
        elif project_name in ["ovarian_cancer_ccca"]:
            # special case
            sample_meta_path = os.path.join(paths["ORIGINAL_DATA_DIR"], f"Data_{sample_id}_Ovarian", "Samples.csv")
            sample_meta_df = pd.read_csv(sample_meta_path)
            sample_meta_df = sample_meta_df.dropna(axis=1, how="all")
            sample_meta_df = sample_meta_df.astype(str)
        else:
            # if not '\t separated (other projects)
            sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])
            sample_meta_df = pd.read_csv(sample_meta_path, usecols=project_cfg["sample_meta_cols"])
        
        # Set sample dicts for projects with matched sample_id to gse_id
        if project_name in ["fetal_gonad", "subtype_evolution", "hgsoc_subtype_define", "hgsoc_tissue_architecture", "gyne_malignant"]:
                        
            sample_dict = dict(zip(sample_meta_df["gse_id"], sample_meta_df["sample_id"]))
            sample_dict_reversed = {v: k for k, v in sample_dict.items()}
        
        # Set marker file path
        if project_name in ["fetal_gonad", "embryos_mixed"]:

            # If dependent on sample_id (i.e. mixed gender)
            markers_file_path = project_cfg["cell_markers_file_path"][sample_id]
            print(markers_file_path)
                
        else: 
            # If not dependent on sample_id (most samples)
            markers_file_path = project_cfg["cell_markers_file_path"]

        try:
            adata_init_path = os.path.join(paths["RAW_DATA_DIR"], f"{sample_id}_init.h5ad")
            adata_filtered_path = os.path.join(paths["RAW_DATA_DIR"], f"{sample_id}_before_qc.h5ad")
            adata_qc_path = os.path.join(paths["RAW_DATA_DIR"], f"{sample_id}_after_qc.h5ad")

            if os.path.exists(adata_qc_path):
                adata_qc = sc.read_h5ad(adata_qc_path)
            else:
                if os.path.exists(adata_filtered_path):
                    adata_filtered = sc.read_h5ad(adata_filtered_path)
                else:
                    if os.path.exists(adata_init_path):
                        adata_init = sc.read_h5ad(adata_init_path)
                    else:
                        
                        if project_name in ["fetal_gonad", "subtype_evolution", "hgsoc_subtype_define", "hgsoc_tissue_architecture", "gyne_malignant"]:
                            gse_id = sample_dict_reversed.get(sample_id)
                            
                            if project_name in ["fetal_gonad"]:

                                adata_init = import_raw_data_fetal_gonad(sample_id, gse_id, paths["ORIGINAL_DATA_DIR"], adata_init_path, sample_meta_df)
                            
                            else:
                            
                                adata_init = import_raw_data_10x(project_name, sample_id, gse_id, paths["ORIGINAL_DATA_DIR"],adata_init_path, sample_meta_df)
                        
                        
                        elif project_name in ["embryos_mixed", "cell_populations"]:    
                            adata_init = import_raw_data_csv(project_name, sample_id, paths["ORIGINAL_DATA_DIR"],adata_init_path, sample_meta_path)
                        
                        elif project_name in ["ovarian_cancer_ccca"]:
                            subproject = sample_id
                            adata_init = import_raw_data_10x_subprojects(subproject, original_data_dir = paths["ORIGINAL_DATA_DIR"], adata_init_path = adata_init_path)
                            
                    adata_filtered = filter_data(sample_id, adata_init, adata_filtered_path, important_genes, CELL_CYCLE_GENES_FILE_PATH, paths["QC_SAVE_DIR"])
                
                adata_qc = run_qc(adata_filtered, adata_qc_path, important_genes)
        except Exception as e:
            print(f"Failed to process {sample_id}: {e}")
            import traceback; traceback.print_exc()
            continue
       

def refine_processed_data(project_name, paths, project_cfg, global_cfg, samples_to_process):

    important_genes = global_cfg["important_genes"]
    cols_to_plot = global_cfg["umap_obs_cols"]
    
    df_cta_gene_expression_list = []

    df_project_cta_gene_expression_save_dir = os.path.join(ALL_RESULTS_DIR, "tables", "cta_genes")
    os.makedirs(df_project_cta_gene_expression_save_dir, exist_ok=True)

    df_project_cta_gene_expression_save_path = os.path.join(
        df_project_cta_gene_expression_save_dir,
        f"{project_name}_cta_gene_expression.csv"
    )

    color_config_path = os.path.join("python_common_scripts", "global_colors.json")

    if os.path.exists(color_config_path):
        data = json.load(open(color_config_path))
        gene_color_map = data.get("gene_color_map", {})
        cell_type_colors = data.get("cell_type_colors", {})
    else:
        gene_color_map = dict(GENE_COLOR_MAP)
        cell_type_colors = dict(CELL_TYPE_COLORS)

    for sample_id in samples_to_process:
        print("Refining sample", sample_id)

        try:

            # Set sample meta path
            # Set sample meta path
            if project_name in ["subtype_evolution", "hgsoc_subtype_define", "hgsoc_tissue_architecture", "gyne_malignant"]:
                # if '\t' separated
                sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])
                sample_meta_df = pd.read_csv(sample_meta_path, usecols=project_cfg["sample_meta_cols"], sep = '\t')
            
            elif project_name in ["ovarian_cancer_ccca"]:
                # special case
                sample_meta_path = os.path.join(paths["ORIGINAL_DATA_DIR"], f"Data_{sample_id}_Ovarian", "Samples.csv")
                sample_meta_df = pd.read_csv(sample_meta_path)
                sample_meta_df = sample_meta_df.dropna(axis=1, how="all")
                sample_meta_df = sample_meta_df.astype(str)
            else:
                # if not '\t separated (other projects)
                sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])
                sample_meta_df = pd.read_csv(sample_meta_path, usecols=project_cfg["sample_meta_cols"])
            
            # Set sample dicts for projects with matched sample_id to gse_id
            if project_name in ["fetal_gonad", "subtype_evolution", "hgsoc_subtype_define", "hgsoc_tissue_architecture", "gyne_malignant"]:
                            
                sample_dict = dict(zip(sample_meta_df["gse_id"], sample_meta_df["sample_id"]))
                sample_dict_reversed = {v: k for k, v in sample_dict.items()}
            
            # Set marker file path
            if project_name in ["fetal_gonad", "embryos_mixed"]:

                # If dependent on sample_id (i.e. mixed gender)
                markers_file_path = project_cfg["cell_markers_file_path"][sample_id]
                print(markers_file_path)
                    
            else: 
                # If not dependent on sample_id (most samples)
                markers_file_path = project_cfg["cell_markers_file_path"]


            adata_qc_path = os.path.join(paths["RAW_DATA_DIR"], f"{sample_id}_after_qc.h5ad")
            adata_annotated_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{sample_id}_annotated_cellmarker.h5ad")
            umap_individual_save_dir = os.path.join(paths["FIGURES_DIR"], "umap", sample_id)

            if os.path.exists(adata_annotated_path):
                adata_annotated = sc.read_h5ad(adata_annotated_path)
            else:
                if not os.path.exists(adata_qc_path):
                    print(f"Missing QC file for {sample_id}, skipping")
                    continue

                adata_qc = sc.read_h5ad(adata_qc_path)

                adata_annotated, gene_color_map, cell_type_colors = cell_annotation_wrapper(
                    sample_id,
                    adata_qc,
                    adata_annotated_path,
                    markers_file_path,
                    gene_color_map,
                    cell_type_colors,
                    genes_list=important_genes,
                    cell_type_label_col="predicted_cell_type",
                    palette="glasbey"
                )

            adata_umap, gene_color_map, cell_type_colors = make_umaps(
                sample_id,
                adata_annotated,
                umap_individual_save_dir,
                gene_color_map,
                cell_type_colors,
                cols_to_plot,
                genes_list=important_genes,
                cell_type_label_col="predicted_cell_type",
                palette="glasbey",
                plot_always=True
            )

            df_cta_gene_expression = cta_genes_expression(
                sample_id,
                adata_umap,
                tables_dir=paths["TABLES_DIR"],
                figures_dir=paths["FIGURES_DIR"],
                cta_genes_file_path=CTA_FAMILY_FILE_PATH,
                top_n=25,
                plot_always=True
            )

            if df_cta_gene_expression is not None:
                df_cta_gene_expression_list.append(df_cta_gene_expression)

            if project_name in ["embryos_mixed"]:
                if "F_" in sample_id:
                    get_ranked_genes(adata = adata_umap, sample_id = sample_id, tables_dir=paths["TABLES_DIR"],
                figures_dir=paths["FIGURES_DIR"], always_rank = True, marker_genes_dict =CUSTOM_MARKER_GENES_DICT)
            elif project_name in ["fetal_gonad"]:
                if "Ovary" in sample_id or "F_Mesonephros" in sample_id:
                    get_ranked_genes(adata = adata_umap, sample_id = sample_id, tables_dir=paths["TABLES_DIR"],
                figures_dir=paths["FIGURES_DIR"], always_rank = True, marker_genes_dict =CUSTOM_MARKER_GENES_DICT)
            else:
                get_ranked_genes(adata = adata_umap, sample_id = sample_id, tables_dir=paths["TABLES_DIR"],
                figures_dir=paths["FIGURES_DIR"],always_rank = True, marker_genes_dict =CUSTOM_MARKER_GENES_DICT)



        except Exception as e:
            print(f"Failed to refine {sample_id}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Combine all samples in the project
    if df_cta_gene_expression_list:
        df_project_cta_gene_expression_df = pd.concat(df_cta_gene_expression_list,ignore_index=True)

        df_project_cta_gene_expression_df.to_csv(df_project_cta_gene_expression_save_path, index=False)
    else:
        print("No CTA gene expression data generated.")

    # Create project-wide CTA dotplot
    if os.path.exists(df_project_cta_gene_expression_save_path):
        all_figures_dir = os.path.join(ALL_RESULTS_DIR, "figures")
        cta_genes_expression_all_samples(project_name,df_project_cta_gene_expression_save_path, all_figures_dir, top_n=25)

    # --- Save colors ---
    os.makedirs(os.path.dirname(color_config_path), exist_ok=True)

    json.dump(
        {
            "gene_color_map": gene_color_map,
            "cell_type_colors": cell_type_colors
        },
        open(color_config_path, "w")
    )