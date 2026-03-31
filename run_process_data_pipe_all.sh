#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --output=logs/run_all/_40G%j.out
#SBATCH --error=errors/run_all/_40G%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=2

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

mkdir -p errors/run_all logs/run_all

cd /data/scottaa/cta_onco_fetal

PROJECTS=(
"cell_populations"
"embryos_mixed"
"fetal_gonad"
"FT_tumors"
#"GSE165897"
#"GSE178101"
"gyne_malignant"
"hgsoc_subtype_define"
"hgsoc_tissue_architecture"
"ovarian_cancer_ccca"
"subtype_evolution"
)

for PROJECT in "${PROJECTS[@]}"; do
    echo "Processing $PROJECT"

    python python_common_scripts/main.py process ${PROJECT} || echo "process failed: $PROJECT"
    python python_common_scripts/main.py refine ${PROJECT} || echo "refine failed: $PROJECT"

done

# To run (example):
# sbatch run_process_data_pipe_all.sh
