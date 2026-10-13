################################################################################
## Run script for data preparation of Chapter 3 (CrossLevelResilience), Baldy
## 
## Written by Hung-wei Lin Oct 2026

rm(list=ls())

# Load packages and source functions -------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

  # Output sub-folders of this site
site_tag <- "Baldy"

# Check the raw macro-environment data (Schofield Pass SNOTEL 737) -------------
  # Show the data coverage of each year, visualized in 'coverage.png'.
  # Flag suspicious data (spike, extreme, flat), visualized in 'series.pdf'.
source(here("code", "Baldy", "CLR_CheckMacroEnv.R"))
