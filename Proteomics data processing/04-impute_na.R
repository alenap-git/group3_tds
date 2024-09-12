# OLINK: Impute missing values (no QC information)

library(BiocManager)
library(impute)

proteins = readRDS("../Data/Proteins_denoised.rds")

imp = impute.knn(as.matrix(proteins),k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed=362436069)
proteins_imp = imp$data

saveRDS(proteins_imp, "../Data/Proteins_na_imputed.rds")