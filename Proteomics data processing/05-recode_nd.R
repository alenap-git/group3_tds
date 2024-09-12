# OLINK: Recoding values below LOD

proteins = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_na_imputed.rds")

dict_qc = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Dictionaries/QC_information_UKB.rds")
dict_sample = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Dictionaries/UKB_sample_information.rds")

## Recode samples with values below LOD ----
X = proteins

all(rownames(X)==dict_sample$eid)

pb = utils::txtProgressBar(style = 3)
LOD_table = NULL
for(k in 1:ncol(X)){
  qc = dict_qc[make.names(dict_qc$Assay)==colnames(X)[k],]
  LOD_list = NULL
  for(i in 1:nrow(X)){
    LOD = qc$LOD[qc$PlateID==dict_sample$PlateID[i]]
    LOD_list = c(LOD_list, LOD) 
  }
  LOD_table = cbind(LOD_table, LOD_list)
  utils::setTxtProgressBar(pb, k/ncol(proteins))
}

X[X < LOD_table] = NA

## Filter assays not detected (<LOD) >=50% of samples ----
ndprop = apply(X, 2, function(x) sum(is.na(x), na.rm = T)/nrow(X))
prop = sort(ndprop, decreasing = TRUE)

pdf("../Figures/Proportion_low_detection.pdf", width = 7, height = 5)
par(mar=c(1,5,2,1), pty = "m")
plot(prop, type = "h", ylab="Proportion of non-detects", xlab="", xaxt = "n", lwd=1, ylim = c(0,1),
     col=ifelse(prop>=0.5, yes="tomato", no="black"), cex.lab=1.5)
abline(h=0.5, lty=2, col="darkred", lwd = 1)
# for (i in 1:length(prop)){
#   axis(side=1, at=i, labels=names(prop)[i], las=2,
#        col=ifelse(prop[i]>=0.5, yes="tomato", no="black"),
#        col.axis=ifelse(prop[i]>=0.5, yes="tomato", no="black"))
# }
dev.off()

proteins_sub = X[,ndprop<0.5]

saveRDS(proteins_sub, "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_nd_filtered.rds")
