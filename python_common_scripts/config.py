# python_scripts_common/config.py

import os
import pandas as pd

BASE_DIR = "/data/scottaa/cta_onco_fetal"

REFERENCE_DIR = os.path.join(BASE_DIR, "reference_data_common")
ALL_RESULTS_DIR = os.path.join(BASE_DIR, "all_results")


ALL_CELL_MARKERS_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cellmarker2_all_sexes_filtered_oncofetal_markers_df.csv")
MALE_CELL_MARKERS_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cellmarker2_male_spec_filtered_oncofetal_markers_df.csv")
FEMALE_CELL_MARKERS_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cellmarker2_female_spec_filtered_oncofetal_markers_df.csv")

ALL_CELL_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cellmarker2_all_sexes_filtered_oncofetal_annotation.csv")
MALE_CELL_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cellmarker2_male_spec_filtered_oncofetal_annotation.csv")
FEMALE_CELL_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cellmarker2_female_spec_filtered_oncofetal_annotation.csv")

ENSEMBL_TO_HGNC_MAP= os.path.join(REFERENCE_DIR,  "ensembl_to_hgnc.csv")
IMPORTANT_GENES = ["CTCF", "CTCFL", "DPEP3"]
GENE_MAP_FILE_PATH = os.path.join(REFERENCE_DIR,"gencode_v45_ensembl_to_hgnc.tsv")
CTA_FAMILY_FILE_PATH = os.path.join(REFERENCE_DIR, "CTA_family.csv")
CELL_CYCLE_GENES_FILE_PATH = os.path.join(REFERENCE_DIR, "cell_cycle_genes.csv")
GENE_DEFAULT_PALETTE = {
    "CTCFL": "#04531c",
    "CTCF": "#4b11d3",
    "DPEP3": "#d31162",
}

GENE_COLOR_MAP = GENE_DEFAULT_PALETTE.copy()
CELL_TYPE_COLORS = {"None": "#d3d3d3"}

# Project Configuration

GLOBAL_CONFIG = {
    "important_genes": IMPORTANT_GENES,
    "global_cell_markers_cols": ['tissue_class', 'tissue_type', 'cancer_type', 'cell_type', 'cell_name',  'Genetype'],
    "umap_obs_cols": ["leiden", "predicted_cell_type"],
    "species": "Homo sapiens", 
    "gene_defaults" : {
        "ctcfl": "#04531c",
        "ctcf":  "#4b11d3",
        "dpep3": "#9B0993",
        "wt1":   "#331704FF",
    },
}

    # Group: fetal

PROJECT_CONFIG = {
    
    "fetal_gonad": {
        
        "sample_col": "sample_id",
        
        "sample_meta_cols": ['sample_id','gse_id', 'Description', 'tissue_source', 'sex', 'trimester', 'number_donors', 'age_gestation'],
        "sample_meta_filename": "fetal_gonad_sample_meta.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "prebuilt_h5",
        "cell_markers_file_path": {
            'G1_1_Ovary' : FEMALE_CELL_MARKERS_FILE_PATH,
            'G2_2_Ovary' :FEMALE_CELL_MARKERS_FILE_PATH,
            'G3_1_Testis' : MALE_CELL_MARKERS_FILE_PATH,
            'G4_2_Testis' :MALE_CELL_MARKERS_FILE_PATH,
            'G_1_2_M_Mesonephros': MALE_CELL_MARKERS_FILE_PATH,
            'G_2_F_Mesonephros':FEMALE_CELL_MARKERS_FILE_PATH,
            'G5_A_2_Ovary' : FEMALE_CELL_MARKERS_FILE_PATH,
            'G5_B_2_Ovary' :FEMALE_CELL_MARKERS_FILE_PATH,
            'G6_A_Mixed': ALL_CELL_MARKERS_FILE_PATH,
            'G6_B_Mixed': ALL_CELL_MARKERS_FILE_PATH,
        },
        
        
        "cell_annotation_file_path":  {
            'G1_1_Ovary' : FEMALE_CELL_ANNOTATION_FILE_PATH,
            'G2_2_Ovary' : FEMALE_CELL_ANNOTATION_FILE_PATH,
            'G3_1_Testis' : MALE_CELL_ANNOTATION_FILE_PATH,
            'G4_2_Testis' :MALE_CELL_ANNOTATION_FILE_PATH,
            'G_1_2_M_Mesonephros': MALE_CELL_ANNOTATION_FILE_PATH,
            'G_2_F_Mesonephros': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'G5_A_2_Ovary' : FEMALE_CELL_ANNOTATION_FILE_PATH,
            'G5_B_2_Ovary' :FEMALE_CELL_ANNOTATION_FILE_PATH,
            'G6_A_Mixed': ALL_CELL_ANNOTATION_FILE_PATH, 
            'G6_B_Mixed': ALL_CELL_ANNOTATION_FILE_PATH,
        },
    },

    "embryos_mixed": {
        
        "sample_col": "sample_id",
        
        "sample_meta_cols": ['sample_id','developmental_stage','sex','cell_type_descript'],
        "sample_meta_filename": "embryos_mixed_concat_sample_meta.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "counts_txt",
        "cell_markers_file_path": {
            'F_10W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_11W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_12W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_14W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_18W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_20W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_23W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_24W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_26W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_5W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_7W': FEMALE_CELL_MARKERS_FILE_PATH,
            'F_8W': FEMALE_CELL_MARKERS_FILE_PATH,
            'M_10W': MALE_CELL_MARKERS_FILE_PATH,
            'M_12W': MALE_CELL_MARKERS_FILE_PATH,
            'M_19W': MALE_CELL_MARKERS_FILE_PATH,
            'M_20W': MALE_CELL_MARKERS_FILE_PATH,
            'M_21W': MALE_CELL_MARKERS_FILE_PATH,
            'M_25W': MALE_CELL_MARKERS_FILE_PATH,
            'M_4W': MALE_CELL_MARKERS_FILE_PATH,
            'M_9W': MALE_CELL_MARKERS_FILE_PATH
        },
        
        "cell_annotation_file_path": {
            'F_10W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_11W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_12W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_14W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_18W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_20W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_23W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_24W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_26W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_5W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_7W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'F_8W': FEMALE_CELL_ANNOTATION_FILE_PATH,
            'M_10W': MALE_CELL_ANNOTATION_FILE_PATH,
            'M_12W': MALE_CELL_ANNOTATION_FILE_PATH,
            'M_19W': MALE_CELL_ANNOTATION_FILE_PATH,
            'M_20W': MALE_CELL_ANNOTATION_FILE_PATH,
            'M_21W': MALE_CELL_ANNOTATION_FILE_PATH,
            'M_25W': MALE_CELL_ANNOTATION_FILE_PATH,
            'M_4W': MALE_CELL_ANNOTATION_FILE_PATH,
            'M_9W': MALE_CELL_ANNOTATION_FILE_PATH
        },
        
    },

    "ovarian_cancer_ccca": {
        "sample_list": ['Geistlinger2020','Izar2020', 'Nath2021','Olalekan2021','Olbrecht2021','Qian2020','Regner2021','Shih2018','Tang-Huau2018','Zhang2019', 'Zhang2022'],
        "sample_col": "sample_id",
        "sample_meta_cols": None,
        "sample_meta_filename": None,
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "10x",
        "cell_markers_file_path": FEMALE_CELL_MARKERS_FILE_PATH,
        "cell_annotation_file_path": FEMALE_CELL_ANNOTATION_FILE_PATH,
    },

}

def set_up_project_config(project_name, base_dir=BASE_DIR):
    # 1. Validation: Use project_name to look up the dictionary key
    if project_name not in PROJECT_CONFIG:
        raise ValueError(f"Project '{project_name}' not found in PROJECT_CONFIG. Available: {list(PROJECT_CONFIG.keys())}")
    project_cfg = PROJECT_CONFIG[project_name]
    global_cfg = GLOBAL_CONFIG
    
    # 2. Path Construction
    project_dir = os.path.join(base_dir, "datasets", project_name)
    paths = {
        "PROJECT_DIR": project_dir,
        "ORIGINAL_DATA_DIR": os.path.join(project_dir, "original_data"),
        "RAW_DATA_DIR": os.path.join(project_dir, "raw_data"),
        "WORKING_ADATA_DIR": os.path.join(project_dir, "working_adata"),
        "FIGURES_DIR": os.path.join(project_dir, "results", "figures"),
        "QC_SAVE_DIR": os.path.join(project_dir, "results", "figures", "qc_plots"),
        "TABLES_DIR": os.path.join(project_dir, "results", "tables")

        
        
    }

    # Auto-create directories
    for path in paths.values():
        os.makedirs(path, exist_ok=True)

    
    # 5. Return everything the pipeline needs
    return {
        "paths": paths,
        "project_cfg": project_cfg,
        "global_cfg": global_cfg,
    }