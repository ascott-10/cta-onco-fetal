#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --output=logs/third_round/_64G%j.out
#SBATCH --error=errors/third_round/_64G%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

mkdir -p errors/third_round logs/third_round
cd /data/scottaa/cta_onco_fetal


PROJECT=$1

python python_common_scripts/main.py process ${PROJECT}
python python_common_scripts/main.py refine ${PROJECT}

# To run (example):
# sbatch run_process_data_pipe.sh embryos_mixed
# sbatch run_process_data_pipe.sh ovarian_cancer_ccca
# sbatch run_process_data_pipe.sh fetal_gonad
# sbatch run_process_data_pipe.sh cell_populations
# sbatch run_process_data_pipe.sh subtype_evolution
# sbatch run_process_data_pipe.sh hgsoc_subtype_define
# sbatch run_process_data_pipe.sh gyne_malignant