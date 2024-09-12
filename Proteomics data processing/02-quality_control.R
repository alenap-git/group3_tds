# OLINK: Filtering / Recoding missing values

library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(stringr)
library(tidyr)
library(renv)

dict_assay = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Dictionaries/Assay_information_UKB.rds")
dict_qc = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Dictionaries/QC_information_UKB.rds")
dict_sample = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Dictionaries/UKB_sample_information.rds")

data = read.csv("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/olink_protein_matrix.csv")
rownames(data) = data$X # eid
proteins = data[,-1]

## Check samples with zero values / missing resources ----
nullprop = apply(proteins, 2, function(x) sum(x==0)/nrow(proteins))
names(nullprop)[which.max(nullprop)]

proteins_recoded = proteins
pb = utils::txtProgressBar(style = 3)
for(k in 1:ncol(proteins)){
  if(sum(proteins[,k]==0)>0){
    null.eid = rownames(proteins)[proteins[,k]==0]
    null.plate = unique(dict_sample$PlateID[which(dict_sample$eid %in% null.eid)])

    tmp = dict_qc[make.names(dict_qc$Assay)==colnames(proteins)[k],]
    tmp2 = tmp[which(tmp$PlateID %in% null.plate),]
    na.plate = tmp2$PlateID[rowSums(is.na(tmp2))>0]
    na.eid = rownames(proteins)[which(dict_sample$PlateID %in% na.plate)]
    proteins_recoded[na.eid,k] = NA
  }
  utils::setTxtProgressBar(pb, k/ncol(proteins))
}

saveRDS(proteins_recoded, "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_na_recoded.rds")

## Check samples with no LOD ----

pb = utils::txtProgressBar(style = 3)
for(k in 1:ncol(proteins)){
  tmp = dict_qc[make.names(dict_qc$Assay)==colnames(proteins)[k],]
  if(sum(is.na(tmp$LOD))>0){
    tmp2 = tmp[is.na(tmp$LOD),]
    na.eid = rownames(proteins)[which(dict_sample$PlateID %in% tmp2$PlateID)]
    naprop = sum(is.na(proteins_recoded[na.eid,k]))/length(na.eid)
    if(naprop!=1){
      print(colnames(proteins)[k])
      print(summary(proteins_recoded[na.eid,k])) 
    }
  }
  utils::setTxtProgressBar(pb, k/ncol(proteins))
}

## Check samples with no Processing Start Date ----

pb = utils::txtProgressBar(style = 3)
for(k in 1:ncol(proteins)){
  tmp = dict_qc[make.names(dict_qc$Assay)==colnames(proteins)[k],]
  if(sum(is.na(tmp$Processing_StartDate))>0){
    tmp2 = tmp[is.na(tmp$Processing_StartDate),]
    na.eid = rownames(proteins)[which(dict_sample$PlateID %in% tmp2$PlateID)]
    naprop = sum(is.na(proteins_recoded[na.eid,k]))/length(na.eid)
    if(naprop!=1){
      print(colnames(proteins)[k])
      print(summary(proteins_recoded[na.eid,k])) 
    }
  }
  utils::setTxtProgressBar(pb, k/ncol(proteins))
}


## Filter assays missing >=50% of samples ----
naprop = apply(proteins_recoded, 2, function(x) sum(is.na(x))/nrow(proteins))
prop = sort(naprop, decreasing = TRUE)
prop = prop[prop>=0.05]
head(prop)

dir.create("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Figures")
pdf("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Figures/Proportion_missing_QC.pdf", width = 7, height = 5)
par(mar=c(5,5,2,1), pty = "m")
plot(prop, type = "h", ylab="Proportion missing QC", xlab="", xaxt = "n", lwd=1, ylim = c(0,1),
     col=ifelse(prop>=0.5, yes="tomato", no="black"), cex.lab=1.5)
abline(h=0.5, lty=2, col="darkred", lwd = 1)
for (i in 1:length(prop)){
  axis(side=1, at=i, labels=names(prop)[i], las=2,
       col=ifelse(prop[i]>=0.5, yes="tomato", no="black"),
       col.axis=ifelse(prop[i]>=0.5, yes="tomato", no="black"))
}
dev.off()

proteins_sub = proteins_recoded[,naprop<0.5]

saveRDS(proteins_sub, "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_na_filtered.rds")

## Check for zero values ----
X = proteins_sub

zeroprop = apply(X, 2, function(x) sum(x==0, na.rm = T)/nrow(X))
summary(zeroprop) # Naturally occurring zero values?

