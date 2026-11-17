################################################################################
## Run script for data preparation of Chapter 3 (CrossLevelResilience), Idaho
## (sagebrush steppe mapped quadrats, Dubois, Idaho, 1923-1957 and 1973)
##
## Written by Hung-wei Lin Oct 2026
##
## Steps up to the crowding kernels are shared by every version of the pipeline and write to the
## 'Idaho' folders; the crowding indices and the Part A input belong to the version set in
## CLR_version_Idaho.R and write to its folders (e.g. 'Idaho(1Env&centreK)').

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

  # Raw data folder of this site, and output sub-folders of the shared steps
site_dir <- "Sagebrush, Dubois, Idaho, 1923-1957"
site_tag <- "Idaho"
cache_dir <- here("data", "data.derived", site_tag, "cache")
  # TRUE: run every step again (e.g. after a change in data.source)
rerun_all <- FALSE

# Check the raw macro-environment data (monthly station climate) ---------------
  # Show the data coverage of each year, visualized in 'coverage.png'.
  # Flag suspicious months (extreme, outlier, annual total mismatch), visualized in 'series.png'.
run_cached(here("code", "Idaho", "CLR_CheckMacroEnv.R"), cache_dir)

# Check the raw population data (mapped quadrats) ------------------------------
  # Show the abundance of every species in every year, visualized in 'abundance.pdf'
  # (size = polygon species, no size = point species, count only = unmapped counts).
  # Compare polygon cover with point plants, visualized in 'pointshare.png'.
run_cached(here("code", "Idaho", "CLR_CheckPop.R"), cache_dir)

# Track individual plants through time (genet IDs) -----------------------------
  # Give every polygon and point a genet ID (trackID) with plantTracker; save the tracked
  # records to 'data/data.derived/Idaho' (about 10 min on the first run, read back later).
run_cached(here("code", "Idaho", "CLR_TrackGenet.R"), cache_dir)

# Check how much the tracking choices change the vital rates -------------------
  # Calibrate the first-year grouping distance against the hand assignment of Chu & Adler (2015),
  # in 'TrackSensitivity_buffGenet.csv'. Re-track the main size species under alternative choices
  # (including the tracked data of Chu & Adler 2015) and compare survival, growth and recruitment
  # estimates in SE units, visualized in 'TrackSensitivity.png' (about 25 min on the first run).
run_cached(here("code", "Idaho", "CLR_TrackSensitivity.R"), cache_dir)

# Build the Part A data (vital-rate model) -------------------------------------
  # One row per genet per year (size, fate at t+1, location, edge, grazing, climate),
  # one row per species, quadrat and annual transition for recruitment,
  # saved to 'data/data.derived/Idaho'; sample sizes per species in 'CLR_PartA_species.csv'.
run_cached(here("code", "Idaho", "CLR_ProcessDataPartA.R"), cache_dir)

# Estimate the crowding kernel and its range -----------------------------------
  # Sum the neighbours of every plant in rings around it (same species, grass, shrub, forb),
  # saved to 'data/data.derived/Idaho' (a few minutes on the first run, read back later).
  # Fit the distance kernels (Teller et al. 2016): one for same species, grass and shrub, one for
  # forb; range R where the crowding effects stop changing, visualized in 'CLR_CrowdKernel_centre.png'
  # (3-4 h on the first run, read back later). Check of one shared kernel against one kernel per
  # growth form (grass, shrub; every penalty) by AIC, in 'CLR_CrowdKernel_split_centre.png'
  # (1-2 h on the first run). Distances from the plant centre; files end in '_centre'.
run_cached(here("code", "Idaho", "CLR_CrowdKernel.R"), cache_dir)

# Estimate the crowding kernel with distances from the plant outline -----------
  # Same as the step above with distances from the outline of the plant (0 inside it), for the
  # "edge" versions; files end in '_edge' (about 2-3 h on the first run). The two kernels compared
  # by AIC in 'CLR_CrowdKernel_origin.csv'.
run_cached(here("code", "Idaho", "CLR_CrowdKernelEdge.R"), cache_dir)

# Build the crowding indices ---------------------------------------------------
  # Crowding of every plant by same species, grass and shrub (main kernel) and by forb (forb
  # kernel), within the range R, in both versions: 'CLR_CrowdIndex_genet_centre.csv' and
  # '_genet_edge.csv', with figures of crowding by plant size ('CLR_CrowdIndex_centre.png',
  # '_edge.png'); quadrat totals for recruitment and the recruitment table of the pooled forb
  # (same in both versions); correlations among the indices in 'CLR_CrowdIndex_cor_centre.csv' (and '_edge').
run_cached(here("code", "Idaho", "CLR_CrowdIndex.R"), cache_dir)

# Build the complete climate record --------------------------------------------
  # Monthly mean temperature and precipitation of the station, 1925-2006 (data set), with missing
  # months filled from the daily records of the same station (GHCN-Daily); seasonal variables of
  # every climate year (July-June): pregrowing season (October-March), growing season (April-June),
  # warmest and coldest month, in 'CLR_Climate_monthly.csv' and 'CLR_Climate_year.csv'.
run_cached(here("code", "Idaho", "CLR_ClimateData.R"), cache_dir)

# Explore the climate variables (PCA) ------------------------------------------
  # PCA of three sets of six seasonal variables (climate year of the transition and the year
  # before), PC1 against PC2 of each set in 'CLR_ClimatePCA.png' (temperature red, precipitation blue);
  # loadings and variance explained in 'CLR_ClimatePCA_loadings.csv'.
run_cached(here("code", "Idaho", "CLR_ClimatePCA.R"), cache_dir)

# Load the workspace of the last shared step -----------------------------------
finish_cached()

# Version of the pipeline ------------------------------------------------------
  # Climate predictor and origin of the crowding distances (CLR_version_Idaho.R); the next step
  # writes to the folders of the version.
source(here("code", "Idaho", "CLR_version_Idaho.R"))
site_tag <- ver_tag
for (f in c("data/data.derived", "result/tables", "result/figs")) dir.create(here(f, site_tag), showWarnings = FALSE)
if (clim_version != "1Env") stop("The climate version '", clim_version, "' is to be written after the climate PCA.")

# Build the input file of the Part A models ------------------------------------
  # Main model (species with at least 500 survival records), forb submodel and other-shrub
  # submodel, with the crowding indices of the version, index tables and scaling constants, in
  # 'CLR_PartA_input.rds'; sample sizes in 'CLR_PartA_input_n.csv'. Skipped when the version's
  # 'CLR_PartA_input.rds' exists (a few minutes).
if (!file.exists(here("data", "data.derived", site_tag, "CLR_PartA_input.rds"))) {
  genet_W <- genet_W_set[[kernel_origin]]
  R_main <- kernel_set[[kernel_origin]]$R_main; R_forb <- kernel_set[[kernel_origin]]$R_forb
  source(here("code", "Idaho", "CLR_InputPartA.R"))
}
