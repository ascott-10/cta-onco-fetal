#!/bin/bash
#SBATCH --job-name=process_cancer
#SBATCH --output=logs/mtab_round/plus_64G%j.out
#SBATCH --error=errors/mtab_round/plus_64G%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

mkdir -p errors/mtab_round logs/mtab_round
cd /data/scottaa/cta_onco_fetal


PROJECT=$1

python python_common_scripts/main.py process ${PROJECT}
#python python_common_scripts/concat_adata.py
#python python_common_scripts/main.py refine ${PROJECT}
#python python_common_scripts/concat_adata.py

# To run (example):
# sbatch run_process_data_pipe.sh mtab_tumors
# sbatch run_process_data_pipe.sh embryos_mixed
# sbatch run_process_data_pipe.sh hgsoc_tumors
# sbatch run_process_data_pipe.sh fetal_gonad