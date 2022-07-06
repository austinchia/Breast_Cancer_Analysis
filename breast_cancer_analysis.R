# About
# This script is used for Breast Cancer Analysis
# This script reads in raw data > preprocesses data > to form a data matrix
# Data matrix is used for further stats analysis

# ============== Clears Memory ======
# clears all objects includes hidden objects
rm(list = ls(all.names = TRUE)) 

# frees up memory and reports the memory usage.
gc() 

# ============== Loads Packages =======
library(readxl)
library(stringr)

# ============== 1. Reads Raw Data ===================

# reads all raw data
bc_raw <- read_excel('Breast_Cancer_Dataset.xlsx', 
                     sheet = '238 spls completed_CBCT_Report', 
                     na = c("", "NA"))

# replaces blank rows with blanks
bc_raw[is.na(bc_raw)] <- ""

# renames columns to show group and subgroup (ie LB_HER2-)
colnames(bc_raw) <- str_replace_all(paste(bc_raw[1,],bc_raw[2,], sep = "_"), '_$','')

# renames first column to accession
colnames(bc_raw)[which(names(bc_raw) == "Group_Subgroup")] <- "Accession"

# uniquefies duplicate column names (LB_1, LB_2)
names(bc_raw) <- make.unique(names(bc_raw), sep="_")

