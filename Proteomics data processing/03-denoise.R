# OLINK: Denoise data using linear mixed models

library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(stringr)
library(tidyr)
library(lme4)



proteins = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_na_filtered.rds")

dict_qc = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Dictionaries/QC_information_UKB.rds")
dict_sample = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Dictionaries/UKB_sample_information.rds")

# Create sample key with plate ID, well ID and processing start date

tech = dict_sample %>% select(PlateID, WellID)
rownames(tech) = dict_sample$eid

# Remove nuisance variance (plate, position and date)

print("Denoising random effects")
pb <- utils::txtProgressBar(style = 3)

denoised = vcov_res = NULL

for (k in 1:ncol(proteins)){
  print(colnames(proteins)[k])
  plate = dict_qc[make.names(dict_qc$Assay)==colnames(proteins)[k],]
  tech_spec = left_join(tech, plate %>% select(PlateID, Processing_StartDate), by = "PlateID")
  rownames(tech_spec) = rownames(tech)
  model1=lmer(proteins[,k]~(1|PlateID)+(1|WellID)+(1|Processing_StartDate), data=tech_spec,
              REML = FALSE,
              control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  vcov=as.data.frame(VarCorr(model1))$vcov
  names(vcov)=as.data.frame(VarCorr(model1))[,1]
  vcov=vcov[-length(vcov)]
  vcov_res = cbind(vcov_res, vcov)
  
  res = residuals(model1)
  names(res) = rownames(proteins)[!is.na(proteins[,k])]
  res = res[rownames(proteins)]
  denoised = cbind(denoised, res)
  
  utils::setTxtProgressBar(pb, k/ncol(proteins))
}
cat("\n")

colnames(denoised) = colnames(proteins)
rownames(denoised) = rownames(proteins)

saveRDS(denoised, "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_denoised.rds")
saveRDS(vcov_res, "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Random_effects_vcov.rds")

