################################################################################
## Run script for the Part A models of Chapter 3 (CrossLevelResilience), Idaho
##
## Written by Hung-wei Lin Oct 2026
##
## Version 1Env&centreK (archived 10 Nov 2026): one climate variable (water-year precipitation) and
## crowding measured from the focal centroid. Outputs in the 'Idaho(1Env&centreK)' folders.
##
## From the full model to model_after. Needs 'CLR_PartA_input.rds' from CLR_runscript_data_Idaho.R.
## Each step is skipped when its output exists (delete the output to run the step again), so the
## whole script can be run again after any step.

rm(list=ls())

# Load packages and source functions -------------------------------------------
library(here)
library(dplyr)
library(coda)     # convergence checks
library(nimble)   # MCMC of the multilevel models (de Valpine et al. 2017)
library(glmmTMB)  # maximum-likelihood fits of the model selection (Brooks et al. 2017)
source(here("code", "CLR_fns_nimble.R"))
source(here("code", "CLR_fns_select.R"))

site_tag <- "Idaho(1Env&centreK)"
fit_path <- function(f) here("data", "data.derived", site_tag, f)

# Full model_before ------------------------------------------------------------
  # Main effects + every interaction term, current precipitation; posterior samples in
  # 'CLR_PartA_stage1_full.rds'.
structure <- "full"
if (!file.exists(fit_path("CLR_PartA_stage1_full.rds"))) source(here("code", "Idaho(1Env&centreK)", "CLR_ModelPartA.R"))

# Submodels: forbs and other shrubs --------------------------------------------
  # Forb (point species pooled) in 'CLR_PartA_forb.rds'; other shrubs (pooled) in 'CLR_PartA_shrub.rds'.
if (!file.exists(fit_path("CLR_PartA_forb.rds"))) source(here("code", "Idaho(1Env&centreK)", "CLR_ModelPartAForb.R"))
if (!file.exists(fit_path("CLR_PartA_shrub.rds"))) source(here("code", "Idaho(1Env&centreK)", "CLR_ModelPartAShrub.R"))

# Plots of the full model ------------------------------------------------------
  # Forest plots in 'PartA_full' ('CLR_PartA_main_effects.png', '_interactions', '_crowding_species',
  # '_climate_grazing_species', '_variance', '_submodels'); Rhat and ESS in 'CLR_PartA_convergence.csv'.
source(here("code", "Idaho(1Env&centreK)", "CLR_PlotPartA.R"))

# Model selection, step 1: climate variable ------------------------------------
  # GLMM fits of the full model with each climate variable (276 fits, about 1-2 hours), shown in
  # 'CLR_PartA_selection_climate.png'.
select_phase <- 1
source(here("code", "Idaho(1Env&centreK)", "CLR_SelectPartA.R"))
if (!phase_done(1)) stop("Some fits of the climate selection failed (errors printed above).")
source(here("code", "Idaho(1Env&centreK)", "CLR_SelectPartA_collect.R"))

# Model selection, step 2: interaction terms (fits on ARC) ---------------------
  # Every subset of interaction terms (66048 fits on ARC, code/arc/CLR_select_array.sh), ranked in
  # 'CLR_PartA_selection_AICcrank.png' and '_BICrank.png'; selected terms in 'CLR_select_selected.csv'.
select_phase <- 2
tasks <- make_tasks(2)
n_missing <- sum(!file.exists(file.path(sel_dir(), tasks$file)))
if (n_missing > 0) stop(sprintf("%d fits of step 2 are missing in 'selection': copy them back from ARC (FileZilla).", n_missing))
rm(tasks, n_missing)
source(here("code", "Idaho(1Env&centreK)", "CLR_SelectPartA_collect.R"))

# Selected model_before (Stage 1) ----------------------------------------------
  # Posterior samples in 'CLR_PartA_stage1_selected.rds'; forest plots in 'PartA_model.selection'.
structure <- "selected"
if (!file.exists(fit_path("CLR_PartA_stage1_selected.rds"))) source(here("code", "Idaho(1Env&centreK)", "CLR_ModelPartA.R"))
source(here("code", "Idaho(1Env&centreK)", "CLR_PlotPartA.R"))

# IPMs and life-history traits (to be written) ---------------------------------

# Model selection, step 3: trait terms of model_after, exploratory (to be written)
# select_phase <- 3
# source(here("code", "Idaho(1Env&centreK)", "CLR_SelectPartA.R"))
# source(here("code", "Idaho(1Env&centreK)", "CLR_SelectPartA_collect.R"))

# model_after (Stage 2; to be written) -----------------------------------------
