################################################################################
## Run script for the Part A models of Chapter 3 (CrossLevelResilience), Idaho
##
## Written by Hung-wei Lin Oct 2026
##
## Needs 'CLR_PartA_input.rds' from CLR_runscript_data_Idaho.R. Each model module has 'test_run'
## in its settings: TRUE for a short check run first. The three models share no parameters.
## A model whose posterior samples already exist is skipped (delete its file to fit it again).

rm(list=ls())

# Load packages and source functions -------------------------------------------
library(here)
library(dplyr)
library(coda)     # convergence checks
library(nimble)   # MCMC of the multilevel models (de Valpine et al. 2017)
source(here("code", "CLR_fns_nimble.R"))

site_tag <- "Idaho"
fit_path <- function(f) here("data", "data.derived", site_tag, f)

# Load or fit the main model, stage 1 (no trait) -------------------------------
  # Joint survival, growth and recruitment model of the main species (baseline for the IPMs
  # that give the life-history trait); posterior samples in 'CLR_PartA_stage1.rds'.
if (!file.exists(fit_path("CLR_PartA_stage1.rds"))) {
  source(here("code", "Idaho", "CLR_ModelPartA.R"))
}

# Load or fit the forb submodel ------------------------------------------------
  # Survival and recruitment of all point species pooled as one forb population;
  # posterior samples in 'CLR_PartA_forb.rds'.
if (!file.exists(fit_path("CLR_PartA_forb.rds"))) {
  source(here("code", "Idaho", "CLR_ModelPartAForb.R"))
}

# Load or fit the other-shrub submodel -----------------------------------------
  # Survival, growth and recruitment of the shrubs below the main-model threshold, pooled;
  # posterior samples in 'CLR_PartA_shrub.rds'.
if (!file.exists(fit_path("CLR_PartA_shrub.rds"))) {
  source(here("code", "Idaho", "CLR_ModelPartAShrub.R"))
}

# Plot the estimates (forest plots) --------------------------------------------
  # Main model: all effects (interaction terms in black; species as faded points) in
  # 'CLR_PartA_main_effects.png', interaction terms alone in 'CLR_PartA_interactions.png', crowding
  # effects per species in 'CLR_PartA_crowding_species.png', climate and grazing per species in
  # 'CLR_PartA_climate_grazing_species.png', variance components in 'CLR_PartA_variance.png'.
  # Submodels (forb, other shrubs) on their own in 'CLR_PartA_submodels.png'.
  # Rhat and effective sample sizes in 'CLR_PartA_convergence.csv'.
source(here("code", "Idaho", "CLR_PlotPartA.R"))
