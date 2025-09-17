# MOCR-DB: Data and Code Repository

This repository provides the core scripts and processed results used in the construction and analysis of **MOCR-DB**.
---
## Repository Structure

- `code/`
  - **app.r** – Shiny application for the MOCR-DB interface
  - **mr_db.r** – Functions for Mendelian randomization (MR) data processing
  - **prepare.r** – Data preparation and preprocessing pipeline

- `data/` (hosted externally due to file size)  
  Processed summary results for MR, LDSC, and SMR analyses, as well as annotation files, are available at:  
  [Google Drive – MOCR-DB data](https://drive.google.com/drive/folders/1JwYUawmOvkdaxxMwNB3gsgqsTTS9SSsP?usp=drive_link)

- `result/`
  - **smr_collect.csv.gz** – Merged SMR result table (required for GO analysis)  
  - **smr_enrich.csv**, **smr_enrich.rdata** – GO enrichment merged output  
  - **data/smr_enrich/** – Per trait GO results
---
## Step 1–4: Pipeline Overview

### Step 1: Heritability & Genetic Correlation (LDSC)
- Scripts:
  - `ldsc_01_munge_h2.sh` – Munge and format input GWAS for LDSC  
  - `ldsc_03_1_ldsc.sh` – Compute single-trait heritability  
  - `ldsc_03_2_ldsc_run.sh` – Pairwise LDSC genetic correlation  
  - `ldsc_4_collect.r` – Merge results across traits  

### Step 2: Mendelian Randomization (MR)
- Scripts:
  - `mr_01_clump.sh` – IV clumping (plink)  
  - `mr_02_pair.r` – Generate trait pairs for MR  
  - `mr_03_mr.r` / `mr_04_mr_run.r` – Perform multiple MR methods  
  - `mr_05_collect.r` – Merge and clean final MR output  

### Step 3: Summary-data-based MR (SMR)
- Scripts:
  - `smr_01_smr.sh` – Run SMR for each QTL-trait combination  
  - `smr_02_collect.r` – Collect & format all SMR results  

### Step 4: GO Enrichment for SMR Results
- Scripts:
  - `smr_00_functions.r` – Internal enrichment functions  
  - `smr_03_enrich.r` – Run enrichment using `smr_collect.csv.gz`  
  - Output stored in `result/smr_enrich/`
---
## Data Sources

- FinnGen: https://www.finngen.fi  
- COVID-19 Host Genetics Initiative: https://www.covid19hg.org/results/r7/  
- UK Biobank and QTL datasets via GCTA Portal:  
  https://yanglab.westlake.edu.cn/software/smr/#DataResource  
---
## Usage

- Scripts in the `code/` directory can be used to process GWAS/QTL data and reproduce key analyses.  
- For GO enrichment, load `smr_00_functions.r` and run `smr_03_enrich.r` using `smr_collect.csv.gz`.  
- Processed datasets (`data/`) provide summary results for integration into the MOCR-DB platform.  
---
## MOCR-DB Platform

The interactive database is publicly available at:  
 https://chenhongwei.net/public/MOCRdb/




