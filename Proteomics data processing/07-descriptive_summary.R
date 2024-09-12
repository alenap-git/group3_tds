# OLINK: Descriptive summary

## Spearman's correlation matrices ----

library(pheatmap)

proteins = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Proteins_nd_imputed.rds")

annot = readRDS("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Data/Panel_annotation.rds")
annotsub = annot[colnames(proteins)]

panel.colour = c("#f4a9be","#61bf1a","#00afac","#9364cc")

mat_col = data.frame(Panel = annotsub)

rownames(mat_col) = colnames(proteins)

mat_colors = list(Panel = panel.colour)

#names(mat_colors$Panel) = levels(annotsub)

panel_levels <- levels(as.factor(mat_col$Panel))

panel_levels

mat_colors <- list(Panel = setNames(panel.colour, panel_levels))

cor = cor(proteins, method = "spearman")
height = 1/144*nrow(cor) + 1
width = height + 2
mybreaks = seq(-1, 1, length = 100)
ifelse(dir.exists("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Figures/"),"",dir.create("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Figures/"))
{png("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Figures/Spearmans_correlation_matrix_hclust.png", width = width, height = height, unit='in', res=500, pointsize = 12)
  pheatmap(cor,
           cellwidth = 0.5, cellheight = 0.5,
           border_color = NA, breaks = mybreaks,
           treeheight_row = 0, treeheight_col = 0,
           show_rownames = FALSE, show_colnames = FALSE,
           annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,
           legend = TRUE, annotation_legend = TRUE,
           annotation_names_row = FALSE, annotation_names_col = FALSE)
  dev.off()
}

{png("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Figures/Spearmans_correlation_matrix.png", width = width, height = height, unit='in', res=500, pointsize = 12)
  pheatmap(cor[order(annotsub),order(annotsub)],
           cellwidth = 0.5, cellheight = 0.5,
           border_color = NA, breaks = mybreaks,
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = FALSE,
           annotation_row = mat_col, annotation_col = mat_col, annotation_colors = mat_colors,
           legend = TRUE, annotation_legend = TRUE,
           annotation_names_row = FALSE, annotation_names_col = FALSE)
  dev.off()
}

tmp = cor
tmp[lower.tri(tmp, diag = T)] = NA
high_corr = cbind(rownames(tmp)[which(abs(tmp) > 0.7, arr.ind = TRUE)[,1]],
                  colnames(tmp)[which(abs(tmp) > 0.7, arr.ind = TRUE)[,2]],
                  tmp[which(abs(tmp) > 0.7, arr.ind = TRUE)])

dir.create("/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Tables")
write.csv(high_corr, "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Tables/Spearmans_correlation_High.csv")

tmp = cor
tmp[lower.tri(tmp, diag = T)] = NA
high_corr2 = cbind(rownames(tmp)[which(abs(tmp) > 0.90, arr.ind = TRUE)[,1]],
                  colnames(tmp)[which(abs(tmp) > 0.90, arr.ind = TRUE)[,2]],
                  tmp[which(abs(tmp) > 0.90, arr.ind = TRUE)])

write.csv(high_corr2, "/rds/general/user/sas123/home/TDS_Project/Data_Preprocessing/Proteomics_data_processing/Tables/Spearmans_correlation_Very_high.csv")
