################################################################################
## Run script for data preparation of Chapter 3 (CrossLevelResilience), Idaho
## (sagebrush steppe mapped quadrats, Dubois, Idaho, 1923-1957 and 1973)
##
## Written by Hung-wei Lin Oct 2026
##
## Version 1Env&centreK (archived 10 Nov 2026): one climate variable (water-year precipitation) and
## crowding measured from the focal centroid. Outputs in the 'Idaho(1Env&centreK)' folders. The saved
## workspaces in its cache hold site_tag "Idaho" (the name before archiving), so to run a step
## again set 'rerun_all <- TRUE'.

rm(list=ls())

# Load packages and source functions -------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)            # shapefiles of the mapped quadrats
library(plantTracker)  # genet tracking (Stears et al. 2022)
  # run_cached(): a module runs once, later runs read its saved workspace until the module
  # script (or a step before it) changes
source(here("code", "CLR_fns_cache.R"))

  # Raw data folder of this site, and output sub-folders of this site
site_dir <- "Sagebrush, Dubois, Idaho, 1923-1957"
site_tag <- "Idaho(1Env&centreK)"
cache_dir <- here("data", "data.derived", site_tag, "cache")
  # TRUE: run every step again (e.g. after a change in data.source)
rerun_all <- FALSE

# Check the raw macro-environment data (monthly station climate) ---------------
  # Show the data coverage of each year, visualized in 'coverage.png'.
  # Flag suspicious months (extreme, outlier, annual total mismatch), visualized in 'series.png'.
run_cached(here("code", "Idaho(1Env&centreK)", "CLR_CheckMacroEnv.R"), cache_dir)

# Check the raw population data (mapped quadrats) ------------------------------
  # Show the abundance of every species in every year, visualized in 'abundance.pdf'
  # (size = polygon species, no size = point species, count only = unmapped counts).
  # Compare polygon cover with point plants, visualized in 'pointshare.png'.
run_cached(here("code", "Idaho(1Env&centreK)", "CLR_CheckPop.R"), cache_dir)

# Track individual plants through time (genet IDs) -----------------------------
  # Give every polygon and point a genet ID (trackID) with plantTracker; save the tracked
  # records to 'data/data.derived/Idaho' (about 10 min on the first run, read back later).
run_cached(here("code", "Idaho(1Env&centreK)", "CLR_TrackGenet.R"), cache_dir)

# Check how much the tracking choices change the vital rates -------------------
  # Calibrate the first-year grouping distance against the hand assignment of Chu & Adler (2015),
  # in 'TrackSensitivity_buffGenet.csv'. Re-track the main size species under alternative choices
  # (including the tracked data of Chu & Adler 2015) and compare survival, growth and recruitment
  # estimates in SE units, visualized in 'TrackSensitivity.png' (about 25 min on the first run).
run_cached(here("code", "Idaho(1Env&centreK)", "CLR_TrackSensitivity.R"), cache_dir)

# Build the Part A data (vital-rate model) -------------------------------------
  # One row per genet per year (size, fate at t+1, location, edge, grazing, climate),
  # one row per species, quadrat and annual transition for recruitment,
  # saved to 'data/data.derived/Idaho'; sample sizes per species in 'CLR_PartA_species.csv'.
run_cached(here("code", "Idaho(1Env&centreK)", "CLR_ProcessDataPartA.R"), cache_dir)

# Estimate the crowding kernel and its range -----------------------------------
  # Sum the neighbours of every plant in rings around it (same species, grass, shrub, forb),
  # saved to 'data/data.derived/Idaho' (a few minutes on the first run, read back later).
  # Fit the distance kernels (Teller et al. 2016): one for same species, grass and shrub, one for
  # forb; range R where the crowding effects stop changing, visualized in 'CLR_CrowdKernel.png'
  # (3-4 h on the first run, read back later).
run_cached(here("code", "Idaho(1Env&centreK)", "CLR_CrowdKernel.R"), cache_dir)

# Build the crowding indices of the Part A model -------------------------------
  # Crowding of every plant by same species, grass and shrub (main kernel) and by forb (forb
  # kernel), within the range R; quadrat totals for recruitment; recruitment table of the pooled
  # forb, saved to 'data/data.derived/Idaho'. Correlations among the indices in 'CLR_CrowdIndex_cor.csv'.
run_cached(here("code", "Idaho(1Env&centreK)", "CLR_CrowdIndex.R"), cache_dir)

# Build the input file of the Part A models ------------------------------------
  # Main model (species with at least 500 survival records), forb submodel and other-shrub
  # submodel, with index tables and scaling constants, in 'CLR_PartA_input.rds';
  # sample sizes in 'CLR_PartA_input_n.csv'.
run_cached(here("code", "Idaho(1Env&centreK)", "CLR_InputPartA.R"), cache_dir)

# Load the workspace of the last step ------------------------------------------
finish_cached()
