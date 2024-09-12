# UKB OLINK Data

UK Biobank (ukb673609) OLINK Proteomics

Please specify path and virtual environment name in job (bash) scripts before submitting jobs.

Data Processing steps (number corresponds to script):
1. Generate data dictionaries and filter based on available data in UKB.
2. Check samples with zero values / missing QC resources (LOD and Processing Start Date); Filter assays missing >=50% of samples: NPM1 and PCOLCE were excluded.
3. Remove nuisance variance (plate, well and processing start date) using linear mixed models.
4. Imputation of missing values using K-nearest neighbour imputation algorithm implemented in R package impute.
5. Assays with >=50% of samples with values below LOD were excluded.
6. Imputation of values below LOD using QRILC (Quantile Regression Imputation of Left-Censored data).
7. Calculate pairwise Spearman's correlation coefficient.
8. Perform Principal Component Analysis (PCA).

Data files required: projects/hda_23-24/live/TDS/Project3/Proteomics_data_processing/Data/
- olink_protein_matrix.csv
- coding143.tsv
- ukb_extracted.rds
  * Number of proteins measured (Field 30900)
  * Plate used for sample run (Field 30901)
  * Well used for sample run (Field 30902)
- olink assay
- olink assay version
- olink assay warning
- olink batch number
- olink limit of detection
- olink panel lot number
- olink processing start date