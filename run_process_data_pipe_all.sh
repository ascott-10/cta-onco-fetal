#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --output=logs/first_round/_40G%j.out
#SBATCH --error=errors/first_round/_40G%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=2

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

cd /data/scottaa/cta_onco_fetal

PROJECTS=("fetal_gonad" "ovarian_cancer_ccca" "embryos_mixed")

for PROJECT in "${PROJECTS[@]}"; do
    echo "Processing $PROJECT"
    python python_common_scripts/main.py process ${PROJECT}
    python python_common_scripts/main.py refine ${PROJECT}
done

# To run (example):
# sbatch run_process_data_pipe.sh
