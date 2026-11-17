################################################################################
## Model selection of Part A (Idaho): one chunk of fits, for an array job on ARC
##
## Written by Hung-wei Lin Oct 2026
##
## Usage (from the project folder): Rscript code/Idaho/CLR_SelectPartA_task.R <phase> <chunk> <chunk_size>
## Runs fits (chunk - 1) * chunk_size + 1 to chunk * chunk_size of the phase (fits already saved are skipped).

args <- as.integer(commandArgs(trailingOnly = TRUE))
suppressPackageStartupMessages({library(here); library(glmmTMB)})
source(here("code", "CLR_fns_select.R"))
site_tag <- "Idaho(1Env&centreK)"
n_fits <- nrow(read.csv(file.path(sel_dir(), sprintf("CLR_select_tasks_phase%d.csv", args[1]))))
for (i in intersect(((args[2] - 1) * args[3] + 1):(args[2] * args[3]), 1:n_fits)) run_task(args[1], i)
