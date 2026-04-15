# python_scripts_common/config.py

import os
import pandas as pd

BASE_DIR = "/data/scottaa/cta_onco_fetal"

REFERENCE_DIR = os.path.join(BASE_DIR, "reference_data_common")
ALL_RESULTS_DIR = os.path.join(BASE_DIR, "all_results")

ALL_CELL_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cm2_updated_all_annotated.csv")
ALL_NONCANCER_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cm2_updated_all_noncancer.csv")
MALE_NONCANCER_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cm2_updated_male_noncancer.csv")
FEMALE_NONCANCER_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cm2_updated_female_noncancer.csv")
FEMALE_CANCER_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "markers", "cm2_updated_adult_female_cancer.csv")

ENSEMBL_TO_HGNC_MAP = os.path.join(REFERENCE_DIR, "ensembl_to_hgnc.csv")
CELL_TYPE_MAPPINGS = os.path.join(REFERENCE_DIR, "JSON", "cell_type_genes_mappings.json")
ADATA_PALETTE_MAP = os.path.join(REFERENCE_DIR, "JSON", "adata_palettes.json")

GENE_MAP_FILE_PATH = os.path.join(REFERENCE_DIR, "gencode_v45_ensembl_to_hgnc.tsv")
CTA_FAMILY_FILE_PATH = os.path.join(REFERENCE_DIR, "CTA_family.csv")
CELL_CYCLE_GENES_FILE_PATH = os.path.join(REFERENCE_DIR, "cell_cycle_genes.csv")

IMPORTANT_GENES = ["CTCF", "CTCFL", "DPEP3"]



GENE_DEFAULT_PALETTE = {
    "CTCFL": "#04531c",
    "CTCF": "#4b11d3",
    "DPEP3": "#d31162",
}

LEIDEN_RES_LIST = [0.5, 1.0,2.0]

GLOBAL_CONFIG = {
    "important_genes": IMPORTANT_GENES,
    "species": "Homo sapiens",
    "leiden_res_list": LEIDEN_RES_LIST,
    "gene_map_file_path": GENE_MAP_FILE_PATH,
    "cell_cycle_genes_path": CELL_CYCLE_GENES_FILE_PATH,
    "cta_genes_path": CTA_FAMILY_FILE_PATH,
    "json_palette_path": ADATA_PALETTE_MAP,
    "json_annotations_path": CELL_TYPE_MAPPINGS,
    
}

PROJECT_CONFIG = {
    
    "fetal_gonad": {
        
        "sample_col": "sample_id",
        "gse_list": ['GSM5506062_G1', 'GSM5506063_G2', 'GSM5506064_G3', 'GSM5506065_G4','GSM5704349_mesonephros', 'GSM6703999_mesonephros_F', 'GSM6704000_G5_A','GSM6704001_G5_B', 'GSM6704002_G6_A', 'GSM6704003_G6_B'],
        "sample_list": ['G1_1_Ovary', 'G2_2_Ovary', 'G3_1_Testis','G4_2_Testis','G_1_2_M_Mesonephros','G_2_F_Mesonephros','G5_A_2_Ovary','G5_B_2_Ovary','G6_A_Mixed','G6_B_Mixed'],
        "group_list": ["Female", "Male", "Mixed"],
        "sample_meta_cols": ['ID_REF','gse_id',	'Description','new_name','tissue_source','sex',	'trimester','number_donors','age_gestation',	'sample_id'],
        "plot": {"obs_cols": ["sample_id", 'age_gestation', 'tissue_source'],
             "obs_titles": ["Sample",  "Age (Gestation Weeks)","Tissue Source",],
             },
        "sample_meta_filename": "fetal_gonad_sample_meta.csv",
        
        "processing_type": "prebuilt_h5",
        
    
    },

    "embryos_mixed": {
        
        "sample_col": "sample_id",
        "sample_list": ['F_3_8W','F_9_12W','F_13_16W','F_17_21W','F_22_26W','M_3_8W','M_9_12W','M_17_21W','M_22_26W'],
        "sample_meta_cols": ['sample_id','developmental_stage','sex','cell_type_descript'],
        "sample_meta_filename": "GSE86146_embryos_mixed_sample_meta.csv",
        "plot": {"obs_cols": ["sample_id", "sex", "stage_label"],
             "obs_titles": ["Sample", "Sex", "Age (Gestation Weeks)"],
             },
        "processing_type": "counts_txt",
        
        
    },
     "mtab_tumors": {
        
        "sample_col": "sample_id",
        #"sample_list": ['38b', '59', '74-1', '79', "mtab_tumors"],
        "sample_list": [ "mtab_tumors"],
        "sample_meta_cols": ['donor_id','sample_id','sex','age','FIGO_grade'],
        "sample_meta_filename": "mtab_tumors_sample_metadata.csv",
        "plot": {"obs_cols": ['sample_id', 'age','FIGO_grade'],
             "obs_titles": ["Sample", "Age", "FIGO_grade"],
             },
        "processing_type": "10x",

        
    },

    "ovarian_cancer_ccca": {
        "sample_list": ['Nath2021','Olalekan2021','Olbrecht2021','Qian2020','Regner2021','Shih2018','Tang-Huau2018','Zhang2019', 'Zhang2022'],
        "sample_col": "sample_id",
        "plot": {"obs_cols": ["sample_id"],
             "obs_titles": ["Sample"],
             },
        "processing_type": "10x",
        "cell_markers_file_path": FEMALE_CANCER_ANNOTATION_FILE_PATH,
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
        "ORIGINAL_ORIGINAL_DATA_DIR": os.path.join(project_dir, "original_original_data"),

        "RAW_DATA_DIR": os.path.join(project_dir, "raw_data"),
        "WORKING_ADATA_DIR": os.path.join(project_dir, "working_adata"),
        
        "FIGURES_DIR": os.path.join(project_dir, "results", "figures"),
        "QC_SAVE_DIR": os.path.join(project_dir, "results", "figures", "qc_plots"),
        
        "TABLES_DIR": os.path.join(project_dir, "results", "tables"),
        

        
        
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