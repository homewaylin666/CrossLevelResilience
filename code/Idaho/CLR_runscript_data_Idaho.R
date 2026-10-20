################################################################################
## Run script for data preparation of Chapter 3 (CrossLevelResilience), Idaho
## (sagebrush steppe mapped quadrats, Dubois, Idaho, 1923-1957 and 1973)
##
## Written by Hung-wei Lin Oct 2026

rm(list=ls())

# Load packages and source functions -------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)            # shapefiles of the mapped quadrats
library(plantTracker)  # genet tracking (Stears et al. 2022)

  # Raw data folder of this site, and output sub-folders of this site
site_dir <- "Sagebrush, Dubois, Idaho, 1923-1957"
site_tag <- "Idaho"

# Check the raw macro-environment data (monthly station climate) ---------------
  # Show the data coverage of each year, visualized in 'coverage.png'.
  # Flag suspicious months (extreme, outlier, annual total mismatch), visualized in 'series.png'.
source(here("code", "Idaho", "CLR_CheckMacroEnv.R"))

# Check the raw population data (mapped quadrats) ------------------------------
  # Show the abundance of every species in every year, visualized in 'abundance.pdf'
  # (size = polygon species, no size = point species, count only = unmapped counts).
  # Compare polygon cover with point plants, visualized in 'pointshare.png'.
source(here("code", "Idaho", "CLR_CheckPop.R"))

# Track individual plants through time (genet IDs) -----------------------------
  # Give every polygon and point a genet ID (trackID) with plantTracker; save the tracked
  # records to 'data/data.derived/Idaho' (about 10 min on the first run, read back later).
source(here("code", "Idaho", "CLR_TrackGenet.R"))

# Check how much the tracking choices change the vital rates -------------------
  # Calibrate the first-year grouping distance against the hand assignment of Chu & Adler (2015),
  # in 'TrackSensitivity_buffGenet.csv'. Re-track the main size species under alternative choices
  # (including the tracked data of Chu & Adler 2015) and compare survival, growth and recruitment
  # estimates in SE units, visualized in 'TrackSensitivity.png' (about 25 min on the first run).
source(here("code", "Idaho", "CLR_TrackSensitivity.R"))

# Build the Part A data (vital-rate model) -------------------------------------
  # One row per genet per year (size, fate at t+1, location, edge, grazing, climate),
  # one row per species, quadrat and annual transition for recruitment,
  # saved to 'data/data.derived/Idaho'; sample sizes per species in 'CLR_PartA_species.csv'.
source(here("code", "Idaho", "CLR_ProcessDataPartA.R"))
