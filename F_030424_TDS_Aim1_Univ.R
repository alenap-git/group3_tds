# Load required packages
library(Matrix)
library(lme4)
library(dplyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(viridis)
suppressPackageStartupMessages(library(VennDiagram))

# Load and preprocess datasets
covars4 <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/covars_2.csv") %>%
  select(-1) %>%
  rename(type = type3) %>%
  mutate(
    type = as.factor(type),
    Sex = as.factor(Sex),
    modified_ethnicity = as.factor(modified_ethnicity),
    Alcohol_intake_recode = as.factor(Alcohol_intake_recode),
    Smoking_status = as.factor(Smoking_status),
    new_IMD_group = as.factor(new_IMD_group),
    comorbidity_category = as.factor(comorbidity_category),
    Controlled_baseline_BP = as.factor(Controlled_baseline_BP)
  ) %>%
  filter(type %in% c("Untreated_ACE/ARB", "ACE/ARBs only")) %>%
  column_to_rownames(var = "X") %>%
  mutate(across(c("Age_at_recruitment", "BMI", "Baseline_systolic_BP"), scale)) %>%
  mutate(type1_ = as.factor(ifelse(type == "Untreated_ACE/ARB", 0, 1)))

# Load protein data
ptns <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/protein_data.csv") %>%
  select(-1) %>%
  column_to_rownames(var = "X") %>%
  scale() %>%
  as.data.frame()

# Subset proteins
ptns4 <- ptns[rownames(covars4), ]

# Define linear model function
foo <- function(X) {
  model0 <- lm(X ~ Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = covars4)
  model1 <- update(model0, . ~ . + type1_)
  res <- c(coef(model1)["type1_1"], pval = anova(model0, model1)$`Pr(>F)`[2])
  return(res)
}

# Apply function and process results
univ4 <- t(apply(ptns4, 2, foo)) %>%
  as.data.frame() %>%
  mutate(pval_bonf = p.adjust(pval, method = "BH"),
         signif = ifelse(pval_bonf < 0.05, "Significant", "Not Significant")) %>%
  rownames_to_column(var = "Protein")

univ4_sign <- filter(univ4, signif == "Significant")
write.table(univ4_sign, "Sign_Hits_ACEs_Unreated-ACEs_FDR.csv", sep = ",", row.names = FALSE)

# Define Volcano plot function
Volcano <- function(results, thr = 0.05, protein_names = NULL) {
  par(mar = c(4.5, 4.5, 1, 1))
  plot(results$type1_1, -log10(results$pval), pch = 19, las = 1, cex = 0.5, 
       xlab = expression(beta), ylab = expression(-log[10](p[value])), 
       col = ifelse(p.adjust(results$pval, method = "BH") < thr, "tomato", "darkgrey"))
  if (!is.null(protein_names)) {
    text(results$type1_1, -log10(results$pval), pos = 3, offset = 0.2, cex = 0.5, labels = protein_names)
  }
  abline(v = 0, lty = 3)
  abline(h = -log10(0.05 / nrow(results)), lty = 2, col = "darkred")
  legend("topleft", col = c("darkred", "tomato", "darkgrey"), lty = c(5, NA, NA), 
         pch = c(NA, 19, 19), cex = 0.7, legend = c("Bonferroni threshold at 0.05", "FDR significant hits", "Not significant"))
}

# Call Volcano plot function
Volcano(univ4)
Ptns_names <- colnames(ptns4)
Volcano(univ4, thr = 0.05 / nrow(univ4), protein_names = Ptns_names)

# Redefine Volcano plot function to exclude max value
Volcano <- function(results, thr = 0.05, protein_names = NULL) {
  par(mar = c(4.5, 4.5, 1, 1))
  max_index <- which.max(-log10(results$pval))
  plot(results$type1_1[-max_index], -log10(results$pval)[-max_index], pch = 19, las = 1, cex = 0.5, 
       xlab = expression(beta), ylab = expression(-log[10](p[value])), 
       col = ifelse(p.adjust(results$pval[-max_index], method = "BH") < thr, "tomato", "darkgrey"))
  if (!is.null(protein_names)) {
    text(results$type1_1[-max_index], -log10(results$pval)[-max_index], pos = 3, offset = 0.2, cex = 0.5, labels = protein_names[-max_index])
  }
  abline(v = 0, lty = 3)
  abline(h = -log10(0.05 / nrow(results)), lty = 2, col = "darkred")
  legend("topleft", col = c("darkred", "tomato", "darkgrey"), lty = c(3, NA, NA), 
         pch = c(NA, 19, 19), cex = 0.7, legend = c("Bonferroni threshold at 0.05", "FDR significant hits", "Not significant"))
}

Volcano(univ4, protein_names = Ptns_names, thr = 0.05 / nrow(univ4))

##################################

# CCB Analysis

# Prepare CCB covariates data
covars5 <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/covars_2.csv") %>%
  select(-1) %>%
  rename(type = type3) %>%
  filter(type %in% c("Untreated_CCB", "CCB only")) %>%
  mutate(
    type = as.factor(type),
    Sex = as.factor(Sex),
    modified_ethnicity = as.factor(modified_ethnicity),
    Alcohol_intake_recode = as.factor(Alcohol_intake_recode),
    Smoking_status = as.factor(Smoking_status),
    new_IMD_group = as.factor(new_IMD_group),
    comorbidity_category = as.factor(comorbidity_category),
    Controlled_baseline_BP = as.factor(Controlled_baseline_BP),
    type1_ = as.factor(ifelse(type == "Untreated_CCB", 0, 1))
  ) %>%
  column_to_rownames(var = "X") %>%
  mutate(across(c("Age_at_recruitment", "BMI", "Baseline_systolic_BP"), scale))

# Subset proteins for CCB
ptns5 <- ptns[rownames(covars5), ]

# Apply linear model function and process results for CCB
univ5 <- t(apply(ptns5, 2, foo)) %>%
  as.data.frame() %>%
  mutate(pval_bonf = p.adjust(pval, method = "BH"),
         signif = ifelse(pval_bonf < 0.05, "Significant", "Not Significant")) %>%
  rownames_to_column(var = "Protein")

univ5_sign <- filter(univ5, signif == "Significant")
write.table(univ5_sign, "Sign_Hits_CCBs_Unreated-CCBs_FDR.csv", sep = ",", row.names = FALSE)

# Call Volcano plot function for CCB
Volcano(univ5)
Ptns_names <- colnames(ptns5)
Volcano(univ5, thr = min(univ5$pval_bonf), protein_names = Ptns_names)
Volcano(univ5, thr = 0.05 / nrow(univ5), protein_names = Ptns_names)

######
# Venn Diagram

# Get significant proteins for ACE and CCB
bonfvl <- rownames(univ4[univ4$pval_bonf < 0.05, ])
bonfl <- rownames(univ5[univ5$pval_bonf < 0.05, ])

# Create Venn diagram
bonf <- list(VL = bonfvl, L = bonfl)
colors <- c("tomato", "forestgreen")
names <- paste0(c("ACEIs/ARB", "CCBs"), " (ptns=", sapply(bonf, length), ")")

venn.diagram(
  bonf,
  filename = "030424_Linear_Univ_Venn_ACE_CCBs_synthctrl_FDR-Labels.png",
  fill = colors,
  category.names = names,
  imagetype = "png",
  cat.pos = c(-90, 90),
  cat.dist = c(0.3, 0.3),
  cat.just = list(c(0, 1.5), c(1.75, 0)),
  cat.cex = 0.6
)
