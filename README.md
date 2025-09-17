# MOCR-DB: Data and Code Repository

This repository provides the core scripts and processed results used in the construction and analysis of **MOCR-DB**.

## Repository Structure

- `code/`
  - **app.r** – Shiny application for the MOCR-DB interface  
  - **mr_db.r** – Functions for Mendelian randomization (MR) data processing  
  - **prepare.r** – Data preparation and preprocessing pipeline  

- `data/` (hosted externally due to file size)  
  Processed summary results for MR, LDSC, and SMR analyses, as well as annotation files, are available at:  
  [Google Drive – MOCR-DB data](https://drive.google.com/drive/folders/1JwYUawmOvkdaxxMwNB3gsgqsTTS9SSsP?usp=drive_link)

## GWAS and QTL Data Sources
- FinnGen: https://www.finngen.fi  
- COVID-19 Host Genetics Initiative: https://www.covid19hg.org/results/r7/  
- UK Biobank and QTL datasets via GCTA Portal: https://yanglab.westlake.edu.cn/software/smr/#DataResource  

## Usage
- Scripts in the `code/` directory can be used to process GWAS/QTL data and reproduce key analyses.  
- Processed datasets (`data/`) provide summary results for integration into the MOCR-DB platform.  

## MOCR-DB Platform
The interactive database is publicly available at:  
https://chenhongwei.net/public/MOCRdb/
