#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --output=logs/first_round/_50G%j.out
#SBATCH --error=errors/first_round/_50G%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

cd /data/scottaa/cta_onco_fetal


# python python_common_scripts/main.py process fetal_gonad
# python python_common_scripts/main.py refine fetal_gonad

python python_common_scripts/main.py process ovarian_cancer_ccca
python python_common_scripts/main.py refine ovarian_cancer_ccca

# To run (example):
sbatch run_process_data_pipe.sh
