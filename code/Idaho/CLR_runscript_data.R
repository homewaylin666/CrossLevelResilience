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
