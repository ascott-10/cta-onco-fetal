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
from codes.import_raw_data import import_raw_data_fetal_gonad, import_raw_data_10x
from codes.filter_qc_data import filter_data, run_qc
from codes.cell_annotation import cell_annotation_wrapper
from codes.visual_gene_expression import make_umaps, cta_genes_expression

def pre_process_project_setup(project_name):
    cfg = set_up_project_config(project_name)
    paths, project_cfg, global_cfg = cfg["paths"], cfg["project_cfg"], cfg["global_cfg"]