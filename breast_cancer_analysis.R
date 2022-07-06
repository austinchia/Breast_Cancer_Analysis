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
library(data.table)
library(tidyverse)

# ============== 1. Reads Raw Data ===================

# reads all raw data
bc_raw <- read_excel('Breast_Cancer_Dataset.xlsx', 
                     sheet = '238 spls completed_CBCT_Report', 
                     na = c("", "NA"))
# ============== 2. Cleans Raw Data =================

# replaces blank rows with blanks
bc_raw[is.na(bc_raw)] <- ""

# renames columns to show group and subgroup (ie LB_HER2-)
colnames(bc_raw) <- str_replace_all(paste(bc_raw[1,],bc_raw[2,], sep = "_"), '_$','')

# renames first column to accession
colnames(bc_raw)[which(names(bc_raw) == "Group_Subgroup")] <- "Accession"

# uniquefies duplicate column names (LB_1, LB_2)
names(bc_raw) <- make.unique(names(bc_raw), sep="_")

# splits multiple accession numbers
bc_raw$Accession <- sapply(strsplit(bc_raw$Accession,";"), `[`, 1)

# ============== 3. Exports Accession Numbers to Upload to Uniprot =====

# exports accession numbers to upload to Uniprot
fwrite(data.frame(bc_raw$Accession), "BC_Accession.csv", sep = ",")

# ============== 4. Combines Uniprot Data To Combined Matrix & Exports Matrix (No Imputation) ======
# reads in Gene Symbol table downloaded from Uniprot
gene_symbol <- fread("BC_Accession_Map.tsv",sep=',')

# splits gene symbol by break
gene_symbol_map <- data.frame(str_split_fixed(gene_symbol$`From	Entry`, '\t',3))

# renames column names
colnames(gene_symbol_map) <- c("Accession", "Entry", "Gene Symbol") 

# splits multiple gene symbols
gene_symbol_split <- data.frame(str_split_fixed(gene_symbol_map$`Gene Symbol`, ' ',4))

# binds gene symbol splits to original gene symbol df
gene_symbol_combined <- cbind(gene_symbol_map, gene_symbol_split)

# selects useful gene symbol columns
gene_symbol_combined <- gene_symbol_combined %>%
  select(c("Accession", "X1")) %>%
  rename("Gene Symbol" = "X1")

# merges gene symbol column to main df
ratio_combined_no_na <- left_join(bc_raw, 
                                  gene_symbol_combined, 
                                  by="Accession") %>%
  relocate(`Gene Symbol`, .after = `Accession`) %>%
  
  # adds number to the end of duplicate gene symbols (ie Sptbn1-2)
  group_by(`Gene Symbol`) %>%
  mutate(`GS_count` = 1:n()) %>%
  mutate(`Gene Symbol` = ifelse(`GS_count` == 1, 
                                `Gene Symbol`, 
                                paste0(`Gene Symbol`, "-", `GS_count`))) %>%
  select(-c(`GS_count`, `Accession`))

# exports combined abundance ratio matrix to csv (no imputation)
fwrite(ratio_combined_no_na, "Output/breast_cancer_combined_GS.csv", sep = ",")

# ============== 5. Imputes Data Using 1/5 of Min. Value & Exports Combined Matrix ======
# removes first 2 rows of group labels
ratio_combined_no_group <- ratio_combined_no_na[-1:-2, ]
rownames(ratio_combined_no_group) <- NULL

# converts dataframe to numeric
ratio_combined_no_group <- ratio_combined_no_group %>%
  mutate_all(function(x) as.numeric(as.character(x)))

# replaces all NAs with 1/5 of minimum positive value
ratio_combined_no_group[2:ncol(ratio_combined_no_group)] <- lapply(ratio_combined_no_group[2:ncol(ratio_combined_no_group)],
                                                                   function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))

# exports combined abundance ratio matrix to csv
fwrite(ratio_combined_no_group, "Output/breast_cancer_combined_GS_imputed.csv", sep = ",")

