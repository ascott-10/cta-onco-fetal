#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --output=logs/run_all/_50G%j.out
#SBATCH --error=errors/run_all/_50G%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

mkdir -p errors/run_all logs/run_all

cd /data/scottaa/cta_onco_fetal

PROJECTS=(
"embryos_mixed"
"fetal_gonad"
)

for PROJECT in "${PROJECTS[@]}"; do
    echo "Processing $PROJECT"

    python python_common_scripts/main.py process ${PROJECT} || echo "process failed: $PROJECT"
    python python_common_scripts/main.py refine ${PROJECT} || echo "refine failed: $PROJECT"

done

# To run (example):
# sbatch run_process_data_pipe_all.sh
