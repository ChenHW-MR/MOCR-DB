# MOCR-DB: Data and Code Repository

This repository contains the key scripts and processed data used for the construction and analysis of MOCR-DB.

## Contents

- `code/`
  - **app.r** – Shiny application script for MOCR-DB interface.  
  - **mr_db.r** – Functions for Mendelian randomization data processing.  
  - **prepare.r** – Data preparation and preprocessing pipeline.  

- `data/`
  - **db_mr.csv** – Processed MR results.  
  - **ldsc_collect.csv** – Processed LDSC results.  
  - **smr_collect.csv** – Processed SMR results.  
  - **tissue_name.csv** – Tissue annotation file.  
  - **link.csv** – Reference links for integrated datasets.  

## Data Sources
- FinnGen: https://www.finngen.fi  
- COVID-19 Host Genetics Initiative: https://www.covid19hg.org/results/r7/  
- UK Biobank and QTL datasets via GCTA Portal: https://yanglab.westlake.edu.cn/software/smr/#DataResource  

## Usage
Scripts in the `code/` folder can be used to process GWAS/QTL data and generate results for visualization in MOCR-DB.  
CSV files in the `data/` folder provide processed summary results for MR, LDSC, and SMR analyses.

## MOCR-DB Platform
The interactive database is available online:  
https://chenhongwei.net/public/MOCRdb/