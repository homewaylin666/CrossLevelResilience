#!/bin/bash
################################################################################
## ARC array job: phase 2 of the model selection of Part A (Idaho), 20 fits per job
##
## Written by Hung-wei Lin Oct 2026
##
## Steps on ARC (from the project folder in $DATA, see the ARC notes in the model runscript):
##   (version of the pipeline as in code/Idaho/CLR_version_Idaho.R; fits go to its 'selection' folder)
##   Rscript -e 'library(here); source(here("code", "Idaho", "CLR_version_Idaho.R")); site_tag <- ver_tag; source(here("code", "CLR_fns_select.R")); nrow(make_tasks(2))'
##   mkdir -p logs
##   sbatch --array=1-10 --export=OFFSET=0 code/arc/CLR_select_array.sh      (test: 10 jobs)
##   sbatch --array=1-1000 --export=OFFSET=0 code/arc/CLR_select_array.sh
##   sbatch --array=1-1000 --export=OFFSET=1000 code/arc/CLR_select_array.sh
##   sbatch --array=1-1000 --export=OFFSET=2000 code/arc/CLR_select_array.sh
##   sbatch --array=1-303 --export=OFFSET=3000 code/arc/CLR_select_array.sh
## 66048 fits = 3303 jobs of 20 fits (about 30-40 min each). Jobs stopped at the time limit can be
## submitted again (saved fits are skipped).

#SBATCH --job-name=clr_select
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/clr_select_%A_%a.out

module load R
cd "$SLURM_SUBMIT_DIR"
Rscript code/Idaho/CLR_SelectPartA_task.R 2 "$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))" 20
