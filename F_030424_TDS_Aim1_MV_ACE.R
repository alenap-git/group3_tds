# Set working directory
setwd("/rds/general/user/sas123/home/TDS_Project/Results_Final_Project/260424_AIM1/Data_Splitting")
getwd()

# 1. Loading libraries
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(fake))
suppressPackageStartupMessages(library(sharp))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(MatchIt))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(epiDisplay))
suppressPackageStartupMessages(library(patchwork))

# 2. Data Loading and Cleaning
covars <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/covars_2.csv")
covars <- covars[, -1, drop = FALSE]
levels(as.factor(covars$type3))
covars$type <- covars$type3
rownames(covars) <- covars[, 1]
covars <- covars[, -1, drop = FALSE]

# 3. Subsetting Data
subset_ACE_2 <- covars$type %in% c("Untreated_ACE/ARB", "ACE/ARBs only")
covars <- covars[subset_ACE_2, ]
covars$type <- as.factor(covars$type)
summary(covars$type)

# 4. Loading Protein Data
proteins <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/protein_data.csv")
proteins <- proteins[, -1, drop = FALSE]
rownames(proteins) <- proteins[, 1]
proteins <- proteins[, -1, drop = FALSE]

# 5. Aligning Protein Data with Covariates
row_names_subset <- rownames(covars)
proteins <- proteins[row_names_subset, ]
covars$Sex <- as.factor(covars$Sex)
covars$modified_ethnicity <- as.factor(covars$modified_ethnicity)
covars$Alcohol_intake_recode <- as.factor(covars$Alcohol_intake_recode)
covars$Smoking_status <- as.factor(covars$Smoking_status)
covars$new_IMD_group <- as.factor(covars$new_IMD_group)
covars$comorbidity_category <- as.factor(covars$comorbidity_category)

# 6. Combining Covariates and Protein Data
combined_data <- cbind(covars, proteins)
summary(combined_data$type)

# 7. Data Splitting into Training and Testing Sets
set.seed(123)
split_index <- createDataPartition(combined_data$type, p = 0.7, list = FALSE)
train_data <- combined_data[split_index, ]
test_data <- combined_data[-split_index, ]

summary(as.factor(train_data$type))
summary(as.factor(test_data$type))

# 8. Data Scaling
common_cols <- intersect(colnames(train_data), colnames(proteins))
additional_cols <- c('Age_at_recruitment', "BMI")
selected_cols <- c(common_cols, additional_cols)

train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)

# 9. Preparing Data for LASSO
train_proteins <- train_data[, intersect(colnames(train_data), colnames(proteins)), drop = FALSE]
train_covars <- train_data[, intersect(colnames(train_data), colnames(covars)), drop = FALSE]
y <- ifelse(train_covars$type == "ACE/ARBs only", 1, 0)

confounders <- c('Age_at_recruitment', "Sex", "comorbidity_category", "Smoking_status", "BMI")
summary(as.factor(train_covars$type))

x <- cbind(train_proteins, train_covars[, confounders])
x <- model.matrix(~ . - 1, data = x)

penalty_factors <- rep(1, ncol(x))
penalty_factors[(ncol(x) - 8):ncol(x)] <- 0

length(penalty_factors) == ncol(x)
v_penalty_factors <- tail(penalty_factors, 15)

# 10. LASSO Cross-Validation
set.seed(1234)
t0 <- Sys.time()
mymodel <- cv.glmnet(x = x, y = y, penalty.factor = penalty_factors, type.measure = "auc", family = "binomial")
t1 <- Sys.time()
print(t1 - t0)
plot(mymodel)

# 11. Selected Variables from LASSO
beta_lasso <- coef(mymodel, s = "lambda.1se")[2:(ncol(proteins) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

# 12. Stability Selection
t0 <- Sys.time()
out <- VariableSelection(xdata = x, ydata = y, verbose = FALSE, penalty.factor = penalty_factors, family = "binomial", n_cat = 3, k = 1000, pi_list = seq(0.5, 0.9, by = 0.01))
t1 <- Sys.time()
print(t1 - t0)
CalibrationPlot(out)
hat_params <- Argmax(out)
print(hat_params)
selprop <- SelectionProportions(out)

par(mar = c(10, 5, 1, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "", ylab = "Selection Proportion", xaxt = "n", col = ifelse(selprop >= hat_params[2], "blue", "grey80"), cex.lab = 1)
abline(h = hat_params[2], lty = 2, col = "red")

threshold_indices <- which(selprop >= hat_params[2])
axis_labels <- names(selprop)[threshold_indices]

for (i in seq_along(threshold_indices)) {
  axis(side = 1, at = threshold_indices[i], labels = axis_labels[i], las = 2, col = "blue", col.axis = "blue", cex.axis = 0.7)
}

selected_protein <- names(selprop)[selprop >= 0.61]
print(selected_protein)
selected_protein <- as.data.frame(selected_protein)

# 13. Logistic Regression Model Training and Testing
train_data$type <- as.factor(train_data$type)
test_data$type <- as.factor(test_data$type)

train_data$type <- relevel(as.factor(ifelse(as.character(train_data$type) == "ACE/ARBs only", yes = 1, no = 0)), ref = "0")
test_data$type <- relevel(as.factor(ifelse(as.character(test_data$type) == "ACE/ARBs only", yes = 1, no = 0)), ref = "0")

summary(train_data$type)

lr.fit <- glm(type ~ ADM+CCL19+CD70+CEP85+CPVL+CRNN+CSF3+DCBLD2+DPEP2+DPP10+DUSP3+FGR+HAVCR1+IL6R+ITGB6+JCHAIN+LAG3+LSP1+MANSC1+MFAP3+MYO9B+MZT1+PPY+PRSS8+REN+S100A16+SFTPD+SLAMF8+SLITRK2+TNFSF11+TYRO3+Age_at_recruitment+Sex+BMI+Smoking_status+comorbidity_category, data = train_data, family = "binomial")

summary(lr.fit)

# Predicting Probabilities and ROC Curve for Training Data
lr.train_pred <- predict(lr.fit, type = "response", newdata = train_data)
train_data <- cbind(train_data, probs = lr.train_pred)

train.pred.obj <- prediction(train_data$probs, train_data$type)
train.perf.obj <- performance(train.pred.obj, "tpr", "fpr")
train.perf.auc <- performance(train.pred.obj, "auc")

plot(train.perf.obj, col = "blue", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.05, labels = paste0("AUC=", round(as.numeric(train.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

train_predicted_class <- factor(ifelse(lr.train_pred > 0.5, "1", "0"), levels = levels(train_data$type))
train_conf_matrix <- confusionMatrix(train_predicted_class, train_data$type)
print(train_conf_matrix)

# Logistic Regression for Testing Data
lr.fit.test <- glm(type ~ ADM+CCL19+CD70+CEP85+CPVL+CRNN+CSF3+DCBLD2+DPEP2+DPP10+DUSP3+FGR+HAVCR1+IL6R+ITGB6+JCHAIN+LAG3+LSP1+MANSC1+MFAP3+MYO9B+MZT1+PPY+PRSS8+REN+S100A16+SFTPD+SLAMF8+SLITRK2+TNFSF11+TYRO3+Age_at_recruitment+Sex+BMI+Smoking_status+comorbidity_category, data = test_data, family = "binomial")

lr.test_pred <- predict(lr.fit.test, type = "response", newdata = test_data)
test_data <- cbind(test_data, probs = lr.test_pred)

test.pred.obj <- prediction(test_data$probs, test_data$type)
test.perf.obj <- performance(test.pred.obj, "tpr", "fpr")
test.perf.auc <- performance(test.pred.obj, "auc")

plot(test.perf.obj, col = "blue", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.05, labels = paste0("AUC=", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

test_predicted_class <- factor(ifelse(lr.test_pred > 0.5, "1", "0"), levels = levels(test_data$type))
test_conf_matrix <- confusionMatrix(test_predicted_class, test_data$type)
print(test_conf_matrix)

# 14. Combined AUC Plot
plot(train.perf.obj, col = "red", main = "", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.06, labels = paste0("Train AUC = ", round(as.numeric(train.perf.auc@y.values), 3)), col = "red", cex = 0.7)
plot(test.perf.obj, col = "blue", add = TRUE)
text(0.9, 0.0, labels = paste0("Test AUC = ", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

# 15. AUC Calculation over 100 Iterations for Test Data
auc_lr_total <- NULL

for (i in 1:100) {
  set.seed(6581 + i)
  split_index <- createDataPartition(combined_data$type, p = 0.7, list = FALSE)
  train_data <- combined_data[split_index, ]
  test_data <- combined_data[-split_index, ]
  
  train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
  test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)
  
  train_data$type <- relevel(as.factor(ifelse(as.character(train_data$type) == "ACE/ARBs only", yes = 1, no = 0)), ref = "0")
  test_data$type <- relevel(as.factor(ifelse(as.character(test_data$type) == "ACE/ARBs only", yes = 1, no = 0)), ref = "0")
  
  lr.fit.test <- glm(type ~ ADM+CCL19+CD70+CEP85+CPVL+CRNN+CSF3+DCBLD2+DPEP2+DPP10+DUSP3+FGR+HAVCR1+IL6R+ITGB6+JCHAIN+LAG3+LSP1+MANSC1+MFAP3+MYO9B+MZT1+PPY+PRSS8+REN+S100A16+SFTPD+SLAMF8+SLITRK2+TNFSF11+TYRO3+Age_at_recruitment+Sex+BMI+Smoking_status+comorbidity_category, data = test_data, family = "binomial")
  
  lr.test_pred <- predict(lr.fit, type = "response", newdata = test_data)
  test_data <- cbind(test_data, probs = lr.test_pred)
  
  test.pred.obj <- prediction(test_data$probs, test_data$type)
  test.perf.obj <- performance(test.pred.obj, "tpr", "fpr")
  test.perf.auc <- performance(test.pred.obj, "auc")
  
  auc_lr <- test.perf.auc@y.values[[1]]
  auc_lr_total <- c(auc_lr_total, auc_lr)
  
  if (i == 1) {
    plot(test.perf.obj, col = "blue", xlab = "1 - Specificity", ylab = "Sensitivity")
  } else {
    plot(test.perf.obj, add = TRUE, col = "blue")
  }
}

abline(a = 0, b = 1, lty = "dashed", col = "gray")
legend("bottomright", legend = paste("Avg Test AUC =", round(mean(auc_lr_total), digits = 3)), col = "blue", lty = 1, cex = 0.8)

# 16. Combined AUC Plot for Training and Testing Data over 100 Iterations
auc_lr_train_total <- NULL
auc_lr_test_total <- NULL

for (i in 1:100) {
  set.seed(6581 + i)
  split_index <- createDataPartition(combined_data$type, p = 0.7, list = FALSE)
  train_data <- combined_data[split_index, ]
  test_data <- combined_data[-split_index, ]
  
  train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
  test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)
  
  train_data$type <- relevel(as.factor(ifelse(as.character(train_data$type) == "ACE/ARBs only", yes = 1, no = 0)), ref = "0")
  test_data$type <- relevel(as.factor(ifelse(as.character(test_data$type) == "ACE/ARBs only", yes = 1, no = 0)), ref = "0")
  
  lr.fit <- glm(type ~ ADM+CCL19+CD70+CEP85+CPVL+CRNN+CSF3+DCBLD2+DPEP2+DPP10+DUSP3+FGR+HAVCR1+IL6R+ITGB6+JCHAIN+LAG3+LSP1+MANSC1+MFAP3+MYO9B+MZT1+PPY+PRSS8+REN+S100A16+SFTPD+SLAMF8+SLITRK2+TNFSF11+TYRO3+Age_at_recruitment+Sex+BMI+Smoking_status+comorbidity_category, data = train_data, family = "binomial")
  
  lr.train_pred <- predict(lr.fit, type = "response", newdata = train_data)
  train_data <- cbind(train_data, probs = lr.train_pred)
  
  train.pred.obj <- prediction(train_data$probs, train_data$type)
  train.perf.obj <- performance(train.pred.obj, "tpr", "fpr")
  train.perf.auc <- performance(train.pred.obj, "auc")
  
  auc_lr_train <- train.perf.auc@y.values[[1]]
  auc_lr_train_total <- c(auc_lr_train_total, auc_lr_train)
  
  if (i == 1) {
    plot(train.perf.obj, col = "red", xlab = "1 - Specificity", ylab = "Sensitivity")
  } else {
    plot(train.perf.obj, add = TRUE, col = "red")
  }
  
  lr.fit.test <- glm(type ~ ADM+CCL19+CD70+CEP85+CPVL+CRNN+CSF3+DCBLD2+DPEP2+DPP10+DUSP3+FGR+HAVCR1+IL6R+ITGB6+JCHAIN+LAG3+LSP1+MANSC1+MFAP3+MYO9B+MZT1+PPY+PRSS8+REN+S100A16+SFTPD+SLAMF8+SLITRK2+TNFSF11+TYRO3+Age_at_recruitment+Sex+BMI+Smoking_status+comorbidity_category, data = test_data, family = "binomial")
  
  lr.test_pred <- predict(lr.fit.test, type = "response", newdata = test_data)
  test_data <- cbind(test_data, probs = lr.test_pred)
  
  test.pred.obj <- prediction(test_data$probs, test_data$type)
  test.perf.obj <- performance(test.pred.obj, "tpr", "fpr")
  test.perf.auc <- performance(test.pred.obj, "auc")
  
  auc_lr_test <- test.perf.auc@y.values[[1]]
  auc_lr_test_total <- c(auc_lr_test_total, auc_lr_test)
  
  if (i == 1) {
    plot(test.perf.obj, col = "blue", xlab = "1 - Specificity", ylab = "Sensitivity")
  } else {
    plot(test.perf.obj, add = TRUE, col = "blue")
  }
  abline(a = 0, b = 1, lty = "dashed", col = "gray")
}

legend(x = "bottom", legend = paste("Avg Training AUC =", round(mean(auc_lr_train_total), digits = 3)), col = "red", lty = 1, box.lty = 0, text.font = 0.1)
legend("bottomright", legend = paste("Avg Test AUC =", round(mean(auc_lr_test_total), digits = 3)), col = "blue", lty = 1, box.lty = 0, text.font = 0.1)

# 17. Forest Plot for Coefficients
extract_coefficients <- function(model_summary) {
  coefficients <- model_summary$coefficients[, 1]
  std_errors <- model_summary$coefficients[, 2]
  odds_ratios <- exp(coefficients)
  ci_lower <- exp(coefficients - 1.96 * std_errors)
  ci_upper <- exp(coefficients + 1.96 * std_errors)
  variables <- rownames(model_summary$coefficients)
  
  df <- data.frame(
    Variable = variables,
    Coefficient = coefficients,
    OR = odds_ratios,
    CI_lower = ci_lower,
    CI_upper = ci_upper
  )
  
  return(df)
}

plot_data1 <- extract_coefficients(summary(lr.fit))
plot_data1$Model <- "Training"
plot_data1$x_max <- max(plot_data1$CI_upper)

plot_data2 <- extract_coefficients(summary(lr.fit.test))
plot_data2$Model <- "Test"
plot_data2$x_max <- max(plot_data2$CI_upper)

combined_plot_data <- rbind(plot_data1, plot_data2)
combined_plot_data$Color <- ifelse(combined_plot_data$CI_lower > 1 | combined_plot_data$CI_upper < 1, "blue", "red")
combined_plot_data$Variable <- factor(combined_plot_data$Variable, levels = rownames(summary(lr.fit.test)$coefficients))
combined_plot_data$Model <- factor(combined_plot_data$Model, levels = c("Training", "Test"))

forest_plot <- ggplot(combined_plot_data, aes(x = OR, y = Variable)) +
  geom_point(aes(color = Color), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0) +
  scale_color_manual(values = c("blue" = "blue", "red" = "red"), guide = "none") +
  labs(title = "", x = "Odds Ratio", y = "Variables") +
  theme_minimal() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  facet_wrap(~Model, scales = "free_x") +
  theme(axis.text.y = element_text(size = 8))

print(forest_plot)

# 18. Additional Forest Plot Customization and Display
new_order <- c("Age_at_recruitment", "SexMale", "BMI", "Smoking_statusNever", "Smoking_statusPrevious", "comorbidity_category1", "comorbidity_category2", "comorbidity_category3+", "ADM", "CCL19", "CD70", "CEP85", "CPVL", "CRNN", "CSF3", "DCBLD2", "DPEP2", "DPP10", "DUSP3", "FGR", "HAVCR1", "IL6R", "ITGB6", "JCHAIN", "LAG3", "LSP1", "MANSC1", "MFAP3", "MYO9B", "MZT1", "PPY", "PRSS8", "REN", "S100A16", "SFTPD", "SLAMF8", "SLITRK2", "TNFSF11", "TYRO3")
new_names <- c("Age", "Male [ref = Female]", "BMI", "Smoking Status = Never [ref = Current]", "Smoking Status = Previous [ref = Current]", "Comorbidity Category 1 [ref = 0]", "Comorbidity Category 2 [ref = 0]", "Comorbidity Category +3 [ref = 0]", "ADM", "CCL19", "CD70", "CEP85", "CPVL", "CRNN", "CSF3", "DCBLD2", "DPEP2", "DPP10", "DUSP3", "FGR", "HAVCR1", "IL6R", "ITGB6", "JCHAIN", "LAG3", "LSP1", "MANSC1", "MFAP3", "MYO9B", "MZT1", "PPY", "PRSS8", "REN", "S100A16", "SFTPD", "SLAMF8", "SLITRK2", "TNFSF11", "TYRO3")

plot_data_list <- list(plot_data1, plot_data2)
modified_data <- reorder_rename_chr(plot_data_list, new_order, new_names)
plot_data1 <- modified_data[[1]]
plot_data2 <- modified_data[[2]]

create_forest_plot <- function(plot_data, model_name) {
  bg_color <- switch(model_name, "Training" = "skyblue", "Test" = "lightpink", "white")
  forest_plot <- ggplot(plot_data, aes(x = OR, y = Variable)) +
    geom_point(shape = 20, size = 4, fill = "black", color = "black") +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0) +
    labs(x = NULL, y = NULL, title = paste(model_name)) +
    theme_minimal() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 9), panel.background = element_rect(fill = bg_color))
  
  if (model_name %in% c("Test")) {
    forest_plot <- forest_plot + theme(axis.text.y = element_blank())
  }
  
  return(forest_plot)
}

plot_model1 <- create_forest_plot(plot_data1, "Training")
plot_model2 <- create_forest_plot(plot_data2, "Test")

combined_plot <- plot_model1 + labs(x = "Odds Ratio") | plot_model2 + labs(x = "Odds Ratio")
print(combined_plot)

# 19. Fold Change Calculation
covars4 <- covars
protein_names <- selected_protein$selected_protein
MV_Sign_ACE <- proteins[, colnames(proteins) %in% protein_names]
MV_Sign_ACE <- rownames_to_column(MV_Sign_ACE, var = "EID")
covars4 <- rownames_to_column(covars4, var = "EID")
MV_Sign_ACE <- merge(MV_Sign_ACE, covars4[, c("EID", "type3")], by = "EID")
names(MV_Sign_ACE)[names(MV_Sign_ACE) == "type3"] <- "Treatment_Status"
summary(as.factor(MV_Sign_ACE$Treatment_Status))

protein_data_long <- melt(MV_Sign_ACE, id.vars = c("EID", "Treatment_Status"))
mean_values <- aggregate(value ~ variable + Treatment_Status, protein_data_long, mean, na.rm = TRUE)

step1_ace_mean <- mean_values[mean_values$Treatment_Status == "ACE/ARBs only", ]
untreated_mean <- mean_values[mean_values$Treatment_Status == "Untreated_ACE/ARB", ]

mean_diff <- merge(step1_ace_mean, untreated_mean, by = "variable", suffixes = c("_Step1_ACE", "_Untreated"), all = TRUE)
mean_diff$Mean_Difference <- mean_diff$value_Step1_ACE - mean_diff$value_Untreated
mean_diff <- mean_diff[order(mean_diff$Mean_Difference, decreasing = TRUE), ]

top_upregulated <- head(mean_diff, 10)
top_downregulated <- tail(mean_diff, 10)

print("Top Upregulated Proteins in Step1 ACE:")
print(top_upregulated)
print("Top Downregulated Proteins in Step1 ACE:")
print(top_downregulated)

threshold <- 0.0
mean_diff$Regulation <- ifelse(mean_diff$Mean_Difference > threshold, "Upregulated", ifelse(mean_diff$Mean_Difference < -threshold, "Downregulated", "Unchanged"))

MV_ACE_Foldchange_Results <- as.data.frame(mean_diff)
write.table(MV_ACE_Foldchange_Results, "MV_ACE_Foldchange_Results.csv", sep = ",", row.names = FALSE)

ggplot(mean_values, aes(x = Treatment_Status, y = variable, fill = value)) +
  geom_tile() +
  labs(x = "Treatment Status", y = "", fill = "NPX Value") +
  scale_fill_viridis_c() +
  scale_x_discrete(labels = c("ACE/ARBs only" = "Step1 ACE/ARBs", "Untreated_ACE/ARB" = "Untreated-ACE/ARBs")) +
  theme_minimal()
