# OLINK data dictionaries

rm(list=ls())
project_path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(project_path)

##install tidyverse dependencies instead 

library(tidyverse)

### Load resources
dat_files_dir = "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data" # Set the directory where your .dat files are located
dat_files = list.files(dat_files_dir, pattern = "\\.dat$", full.names = TRUE) # List all .dat files in the directory
dat_list = list() # Initialize an empty list to store data frames

codings = read.delim("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/coding143.tsv")
extracted = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/ukb_extracted.rds")

# Loop through each .dat file, read it, and store it in the list
for (dat_file in dat_files) {
  # Assuming tab-separated data; adjust the separator if necessary
  dat = read.table(dat_file, header = TRUE, sep = "\t")
  
  # Store the data frame in the list
  dat_list[[basename(dat_file)]] = dat
}

head(dat_list[[1]])
head(dat_list[[2]])
head(dat_list[[3]])
head(dat_list[[4]])
head(dat_list[[5]])
head(dat_list[[6]])
head(dat_list[[7]])

table(dat_list[[2]]$Assay_Warning) # Assays that do not pass the QC are not included

### Create dictionaries
assay = merge(dat_list[[3]],dat_list[[1]])
dict_assay = assay %>% select(Assay, UniProt, Panel) %>% cbind(codings)

number = merge(dat_list[[3]],dat_list[[6]], by = "Panel") %>%
  merge(dat_list[[4]], by = "Batch")
date = merge(number, dat_list[[7]], by = c("PlateID","Panel"), all = TRUE)
summary(date)

table(date$Batch[is.na(date$Processing_StartDate)])
table(date$PlateID[is.na(date$Processing_StartDate)])

lod = dat_list[[5]]
dict_qc  = merge(lod[lod$Instance==0,], date, by = c("PlateID","Assay"), all = TRUE)
dict_qc$PlateID = paste0(0,dict_qc$PlateID)

dict_sample = extracted[which(extracted$`30900-0.0`>0),]
colnames(dict_sample) = c("eid","Number of proteins measured","PlateID","WellID")

dir.create("Dictionaries")
saveRDS(dict_assay, "Dictionaries/Assay_information.rds")
saveRDS(dict_qc, "Dictionaries/QC_information.rds")
saveRDS(dict_sample, "Dictionaries/UKB_sample_information.rds")

## Filter/Amend dictionary based on data set----

dict_qc = dict_qc[which(dict_qc$Panel_Lot_Nr %in% c("B04411","B04412","B04413","B04414")),]
dict_qc = dict_qc[dict_qc$PlateID %in% dict_sample$PlateID,]
dict_qc  = dict_qc %>% mutate_if(is.factor,droplevels)

saveRDS(dict_qc, "Dictionaries/QC_information_UKB.rds")

dict_assay = dict_assay[order(dict_assay$Panel),]
keep = unique(dict_qc$Assay)
dict_assay  = dict_assay %>% filter(Assay %in% keep) %>% mutate_if(is.factor,droplevels)

dict_assay$column = make.names(dict_assay$Assay)

annot = dict_assay$Panel
names(annot) = dict_assay$column

saveRDS(dict_assay, "Dictionaries/Assay_information_UKB.rds")
saveRDS(annot, "Data/Panel_annotation.rds")

