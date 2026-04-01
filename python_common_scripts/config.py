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

CUSTOM_MARKER_GENES_DICT =  {
    "PGC": ["POU5F1", "TCL1A", "PDPN", "MEG3", "GPHA2"],
    "Pre-meiotic": ["DAZL", "ZGLP1", "DDX4", "SMC1B", "TEX30", "CKS2"],
    "Meiotic": ["SYCP3", "SYCP2", "MEIOB", "SPATA22", "ZCWPW1", "MAEL"],
    "Oocyte": ["ZP3", "TDRD1", "BRDT", "RBM46"],
    "Endothelial": ["CDH5", "PECAM1", "VWF", "KDR", "FLT1"],
    "Fibroblast": ["COL1A1", "COL3A1", "DCN", "PDGFRA"]
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
        "sample_list": ['G1_1_Ovary', 'G2_2_Ovary', 'G3_1_Testis','G4_2_Testis','G_1_2_M_Mesonephros','G_2_F_Mesonephros','G5_A_2_Ovary','G5_B_2_Ovary','G6_A_Mixed','G6_B_Mixed'],
        "sample_meta_cols": ['sample_id','gse_id', 'Description', 'tissue_source', 'sex', 'trimester', 'number_donors', 'age_gestation'],
        "sample_meta_filename": "fetal_gonad_sample_meta.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "prebuilt_h5",
        "cell_markers_file_path": {
            'G1_1_Ovary' : FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'G2_2_Ovary' : FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'G3_1_Testis' : MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'G4_2_Testis' :MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'G_1_2_M_Mesonephros': MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'G_2_F_Mesonephros':FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'G5_A_2_Ovary' : FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'G5_B_2_Ovary' :FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'G6_A_Mixed': ALL_NONCANCER_ANNOTATION_FILE_PATH,
            'G6_B_Mixed': ALL_NONCANCER_ANNOTATION_FILE_PATH,
        },
    
    },

    "embryos_mixed": {
        
        "sample_col": "sample_id",
        "sample_list": ['F_3_8W','F_9_12W','F_13_16W','F_17_21W','F_22_26W','M_10W','M_3_8W','M_9_12W','M_20W','M_21W','M_17_21W','M_22_26W'],
        "sample_meta_cols": ['sample_id','developmental_stage','sex','cell_type_descript'],
        "sample_meta_filename": "embryos_mixed_concat_sample_meta.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "counts_txt",
        "cell_markers_file_path": {
            'F_3_8W': FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'F_9_12W': FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'F_13_16W': FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'F_17_21W': FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'F_22_26W': FEMALE_NONCANCER_ANNOTATION_FILE_PATH,
            'M_10W': MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'M_3_8W': MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'M_9_12W': MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'M_20W': MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'M_21W':MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'M_17_21W': MALE_NONCANCER_ANNOTATION_FILE_PATH,
            'M_22_26W': MALE_NONCANCER_ANNOTATION_FILE_PATH,
    
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
        "cell_markers_file_path": FEMALE_CANCER_ANNOTATION_FILE_PATH,
    },

     "cell_populations": {
        "sample_list":  ['349', '553', '565', '568', '580', '589', '600', '618', '626'],

        "sample_col": "subject_id",
        "sample_meta_cols": ["sample_id","subject_id","gse_id","lesion_diagnosis","lesion_site","cancer_type"],
        "sample_meta_filename": "cell_populations_sample_meta.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "counts_txt",
        "cell_markers_file_path": FEMALE_CANCER_ANNOTATION_FILE_PATH,
    },

    "subtype_evolution": {
        "sample_list":  ['T59', 'T76', 'T77', 'T89', 'T90'],

        "sample_col": "sample_id",
        "sample_meta_cols": ["sample_id","gse_id","donor_id","sample_tissue"],
        "sample_meta_filename": "subtype_evolution_metadata.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "10x",
        "cell_markers_file_path": FEMALE_CANCER_ANNOTATION_FILE_PATH,
    },

    "hgsoc_subtype_define": {
        "sample_list":  ['16030X2_HJVMLDMXX', '16030X3_HJTWLDMXX', '16030X4_HJTWLDMXX'],
        "sample_col": "sample_id",
        "sample_meta_cols": ["sample_id","gse_id","sample_tissue","sample_tissue_description"],
        "sample_meta_filename": "hgsoc_subtype_define_sample_meta.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "10x",
        "cell_markers_file_path": FEMALE_CANCER_ANNOTATION_FILE_PATH,
    },

    "hgsoc_tissue_architecture": {
        "sample_list":  ['Norm1','Norm2','Norm3','Norm4','Norm5','Cancer1','Cancer2','Cancer3','Cancer4','Cancer5','Cancer6','Cancer7'],
        "sample_col": "sample_id",
        "sample_meta_cols": ["sample_id","gse_id","sample_tissue","sample_tissue_description", "tumor_stage","sample_long_id"],
        "sample_meta_filename": "hgsoc_tissue_architecture_sample_meta.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
        "processing_type": "10x",
        "cell_markers_file_path": FEMALE_CANCER_ANNOTATION_FILE_PATH,
    },

    "gyne_malignant":
    {
        "sample_list":  ['3533EL','3571DL','36186L','36639L','366C5L','37EACL','38FE7L','3BAE2L','3CCF1L','3E4D1L','3E5CFL'],
        "sample_col": "sample_id",
        "sample_meta_cols": ["sample_id","gse_id","sample_tissue","sample_tissue_description", "donor_ethnicity", "tumor_site","tumor_histology","metastatic_status"],
        "sample_meta_filename": "gynecological_malignancies_sample_meta.csv",
        "cell_meta_cols": None,
        "cell_meta_filename": None,
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