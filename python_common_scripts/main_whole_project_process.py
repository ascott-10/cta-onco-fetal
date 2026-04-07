# python_scripts_common/main_process_data.py

def stop():
        umap_individual_save_dir = os.path.join(paths["FIGURES_DIR"], "umap", sample_id)
        
        
        df_cta_gene_expression_list = []

        df_project_cta_gene_expression_save_dir = os.path.join(ALL_RESULTS_DIR, "tables", "cta_genes")
        os.makedirs(df_project_cta_gene_expression_save_dir, exist_ok=True)

        df_project_cta_gene_expression_save_path = os.path.join( df_project_cta_gene_expression_save_dir, f"{project_name}_cta_gene_expression.csv")

        color_config_path = os.path.join("python_common_scripts", "global_colors.json")
        if os.path.exists(color_config_path):
        data = json.load(open(color_config_path))
        gene_color_map = data.get("gene_color_map", {})
        cell_type_colors = data.get("cell_type_colors", {})
    else:
        gene_color_map = dict(GENE_COLOR_MAP)
        cell_type_colors = dict(CELL_TYPE_COLORS)

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
                
                if project_name in ["ovarian_cancer_ccca"]:
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
                if project_name in ["fetal_gonad"]:
                                
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
                    plot_always=False
                )

                df_cta_gene_expression = cta_genes_expression(
                    sample_id,
                    adata_umap,
                    tables_dir=paths["TABLES_DIR"],
                    figures_dir=paths["FIGURES_DIR"],
                    cta_genes_file_path=CTA_FAMILY_FILE_PATH,
                    top_n=25,
                    plot_always=False
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
        
    def process_whole_project(project_name, paths, global_cfg):

        cluster_celltype_mapping_json = CELL_TYPE_MAPPINGS
        marker_genes_list = []
        for k,v in CUSTOM_MARKER_GENES_DICT.items():
            for genes in v:
                marker_genes_list.append(genes)

        marker_genes_list = list(set(marker_genes_list))

        
        all_figures_dir = os.path.join(ALL_RESULTS_DIR, "figures")
        all_tables_dir = os.path.join(ALL_RESULTS_DIR, "tables")

        
        adata_concat_female_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_female_processed.h5ad")
        adata_concat_male_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_male_processed.h5ad")

        color_config_path = os.path.join("python_common_scripts", "global_colors.json")

        if os.path.exists(color_config_path):
            data = json.load(open(color_config_path))
            gene_color_map = data.get("gene_color_map", {})
            cell_type_colors = data.get("cell_type_colors", {})
        else:
            gene_color_map = dict(GENE_COLOR_MAP)
            cell_type_colors = dict(CELL_TYPE_COLORS)

        df_cta_gene_expression_list = []
        df_project_cta_gene_expression_save_dir = os.path.join(all_tables_dir, "cta_genes")
        os.makedirs(df_project_cta_gene_expression_save_dir, exist_ok=True)
        df_project_cta_gene_expression_save_path = os.path.join( df_project_cta_gene_expression_save_dir, f"{project_name}_cta_gene_expression.csv")


        try:
        
            umap_dir = os.path.join(all_figures_dir,"marker_genes","umap")
            if os.path.exists(adata_concat_female_path):
                adata_concat_female = sc.read_h5ad(adata_concat_female_path)
                dataset_name = f"{project_name}_female"
                umap_dir = os.path.join(all_figures_dir,"marker_genes","umap", dataset_name)
                #get_ranked_genes(adata = adata_concat_female, sample_id = project_name, tables_dir = all_tables_dir, figures_dir=all_figures_dir,always_rank = True, marker_genes_dict =CUSTOM_MARKER_GENES_DICT, save_name="female")
                #adata_concat_female = plot_marker_genes_leiden(adata = adata_concat_female, adata_path = adata_concat_female_path, cell_type_colors = cell_type_colors, gene_color_map = gene_color_map, leiden_res=0.50)
                #adata_concat_female, gene_color_map, cell_type_colors = make_umaps(sample_id = dataset_name, adata = adata_concat_female, umap_dir = umap_dir, gene_color_map = gene_color_map, cell_type_colors = cell_type_colors, cols_to_plot = ["leiden_res_0.50", "sample_id", "predicted_cell_type"], genes_list = marker_genes_list, cell_type_label_col="predicted_cell_type", palette="glasbey", plot_always=True)
                adata_concat_female = rank_genes_for_publication(adata_female_path = adata_concat_female_path, project_name = "embryos_mixed", final_leiden_res = "2.0", cluster_celltype_mapping_json = cluster_celltype_mapping_json, cta_genes_file_path = CTA_FAMILY_FILE_PATH, figures_dir = all_figures_dir)
            #if os.path.exists(adata_concat_male_path):
                #dataset_name = f"{project_name}_male"
                #umap_dir = os.path.join(all_figures_dir,"marker_genes","umap", dataset_name)
                #adata_concat_male = sc.read_h5ad(adata_concat_male_path)
                #get_ranked_genes(adata = adata_concat_male, sample_id = project_name, tables_dir = all_tables_dir, figures_dir=all_figures_dir,always_rank = True, marker_genes_dict =CUSTOM_MARKER_GENES_DICT,save_name="male")
                #adata_concat_male = plot_marker_genes_leiden(adata = adata_concat_male, adata_path = adata_concat_male_path, cell_type_colors = cell_type_colors, gene_color_map = gene_color_map, leiden_res=0.50)
                #adata_concat_male, gene_color_map, cell_type_colors = make_umaps(sample_id = dataset_name, adata = adata_concat_male, umap_dir = umap_dir, gene_color_map = gene_color_map, cell_type_colors = cell_type_colors, cols_to_plot = ["leiden_res_0.50", "sample_id", "predicted_cell_type"], genes_list = marker_genes_list, cell_type_label_col="predicted_cell_type", palette="glasbey", plot_always=True)
            
            if df_cta_gene_expression_list:
                df_project_cta_gene_expression_df = pd.concat(df_cta_gene_expression_list,ignore_index=True)
                df_project_cta_gene_expression_df.to_csv(df_project_cta_gene_expression_save_path, index=False)
            else:
                print("No CTA gene expression data generated.")
            
            if os.path.exists(df_project_cta_gene_expression_save_path):
                
                cta_genes_expression_all_samples(project_name,df_project_cta_gene_expression_save_path, all_figures_dir, top_n=25)
        
        except Exception as e:
                print(f"Failed to continue processing project {project_name}: {e}")
                import traceback
                traceback.print_exc()

        # --- Save colors ---
        os.makedirs(os.path.dirname(color_config_path), exist_ok=True)

        json.dump(
            {
                "gene_color_map": gene_color_map,
                "cell_type_colors": cell_type_colors
            },
            open(color_config_path, "w")
        )


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

    for sample_id in samples_to_process:
        print("Processing sample", sample_id)

        if project_name in ["ovarian_cancer_ccca"]:
            # special case
            sample_meta_path = os.path.join(paths["ORIGINAL_DATA_DIR"], f"Data_{sample_id}_Ovarian", "Samples.csv")
            sample_meta_df = pd.read_csv(sample_meta_path)
            sample_meta_df = sample_meta_df.dropna(axis=1, how="all")
            sample_meta_df = sample_meta_df.astype(str)
        else:
            # if not '\t separated (other projects)
            sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])
            sample_meta_df = pd.read_csv(sample_meta_path, usecols=project_cfg["sample_meta_cols"])
        
        
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
            adata_filt_path = os.path.join(paths["RAW_DATA_DIR"], f"{sample_id}_before_qc.h5ad")
            adata_qc_path = os.path.join(paths["RAW_DATA_DIR"], f"{sample_id}_after_qc.h5ad")

            if os.path.exists(adata_qc_path):
                adata_qc = sc.read_h5ad(adata_qc_path)
            else:
                if os.path.exists(adata_filt_path):
                    adata_filtered = sc.read_h5ad(adata_filt_path)
                else:
                    if os.path.exists(adata_init_path):
                        adata_init = sc.read_h5ad(adata_init_path)
                    else:
                        
                        if project_name in ["fetal_gonad"]:
                            sample_dict = dict(zip(sample_meta_df["gse_id"], sample_meta_df["sample_id"]))
                            sample_dict_reversed = {v: k for k, v in sample_dict.items()}
                            gse_id = sample_dict_reversed.get(sample_id)
                            adata_init = import_raw_data_fetal_gonad(sample_id, gse_id, paths["ORIGINAL_DATA_DIR"], adata_init_path, sample_meta_df)
                        
                        elif project_name in ["embryos_mixed"]:
                            original_original_data_dir = os.path.join(paths["ORIGINAL_ORIGINAL_DATA_DIR"])
                            original_data_dir = os.path.join(paths["ORIGINAL_DATA_DIR"])
                            sample_meta_path = os.path.join(paths["PROJECT_DIR"], project_cfg["sample_meta_filename"])
                            adata_init = import_raw_data_embryos_mixed(project_name, original_original_data_dir, original_data_dir, adata_init_path, sample_meta_path)
                        
    
                        elif project_name in ["ovarian_cancer_ccca"]:
                            subproject = sample_id
                            markers_file_path = project_cfg["cell_markers_file_path"]
                            adata_init = import_raw_data_10x_subprojects(subproject, original_data_dir = paths["ORIGINAL_DATA_DIR"], adata_init_path = adata_init_path)
                    if project_name in ["embryos_mixed"]:
                        adata_filtered = filter_data_embryos_mixed(project_name, adata_init, adata_filt_path, important_genes, cell_cycle_genes_file_path, paths["QC_SAVE_DIR"])   
                    else:

                        adata_filtered = filter_data(sample_id, adata_init, adata_filt_path, important_genes, cell_cycle_genes_file_path, paths["QC_SAVE_DIR"])
                    
                    
                if project_name in ["embryos_mixed"]:
                    adata_qc = run_qc_embryos_mixed(project_name, adata_filtered, adata_qc_path, important_genes, paths["QC_SAVE_DIR"])
                else:
                    adata_qc = run_qc(adata_filtered, adata_qc_path, important_genes)
        
        except Exception as e:
            print(f"Failed to process {sample_id}: {e}")
            import traceback; traceback.print_exc()
            continue
       

def refine_processed_data(project_name, paths, project_cfg, global_cfg, samples_to_process):

    important_genes = global_cfg["important_genes"]
    # Variable setup
    final_leiden_res = global_cfg["final_leiden_res"]
  
    subproject_name = f'{project_name}_{global_cfg["sex"]}'
    
    top_markers_path = os.path.join(paths["TABLES_DIR"], "deg_analysis",global_cfg["top_markers_filename"])
    cta_genes_path = global_cfg["cta_genes_path"]
    json_annotations_path =  global_cfg["json_annotations_path"]

    
    


    # Add cell annotations to individual projects
    if project_name in ["embryos_mixed"]:

        print("Refining project", project_name)

        cols_to_plot = ["sample_id", "sex", "stage_label"]
        cols_to_plot_titles = ["Sample", "Sex", "Age (Gestation Weeks)"]

        adata_qc_path = os.path.join(paths["RAW_DATA_DIR"], f"{project_name}_mixed_after_hvg.h5ad")
        adata_annotated_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_mixed_processed.h5ad")
        adata_male_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_male_processed.h5ad")
        adata_female_path = os.path.join(paths["WORKING_ADATA_DIR"], f"{project_name}_female_processed.h5ad")

        try:
            if os.path.exists(adata_annotated_path):
                adata_annotated = sc.read_h5ad(adata_annotated_path)
            
                adata_female, adata_male = adata_split_annotate(subproject_name, final_leiden_res, adata_female_path, adata_male_path, adata_annotated_path, json_annotations_path)
                adata_female = annotate_obs_columns(subproject_name, adata_female, cols_to_plot, cols_to_plot_titles, adata_female_path, palette = "tab10")
                markers_leiden_embryos_mixed(adata_female, subproject_name, top_markers_path, final_leiden_res, paths["FIGURES_DIR"])
            else:
                print("No QC adata available")

        except Exception as e:
            print(f"Failed to refine {project_name}: {e}")
            import traceback
            traceback.print_exc()
            #continue
    