# OLINK: Impute values below LOD


library(BiocManager)
library(imputeLCMD)

proteins = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_nd_filtered.rds")

imp = impute.QRILC(proteins, tune.sigma = 1)
proteins_imp = imp[[1]]

saveRDS(proteins_imp, "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_nd_imputed.rds")
