# Setting the Working Directory
setwd("/rds/general/user/sas123/home/TDS_Project/Results_Final_Project/260424_AIM1/Data_Splitting")
getwd()

# 1. Loading Required Libraries
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
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))

# 2. Data Loading and Preprocessing
## Loading Covariate Data
covars <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/covars_2.csv")
covars <- covars[, -1, drop = FALSE]

## Renaming and Filtering Columns
covars$type <- covars$type3
rownames(covars) <- covars[, 1]
covars <- covars[, -1, drop = FALSE]

## Filtering Specific Treatment Types
subset_CCB_2 <- covars$type %in% c("Untreated_CCB", "CCB only")
covars <- covars[subset_CCB_2, ]
covars$type <- as.factor(covars$type)

## Loading Protein Data
proteins <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/protein_data.csv")
proteins <- proteins[, -1, drop = FALSE]
rownames(proteins) <- proteins[, 1]
proteins <- proteins[, -1, drop = FALSE]

## Filtering Protein Data Based on Covariate Row Names
row_names_subset <- rownames(covars)
proteins <- proteins[row_names_subset, ]

# 3. Data Integration
## Combining Covariates and Protein Data
covars$Sex <- as.factor(covars$Sex)
covars$modified_ethnicity <- as.factor(covars$modified_ethnicity)
covars$Alcohol_intake_recode <- as.factor(covars$Alcohol_intake_recode)
covars$Smoking_status <- as.factor(covars$Smoking_status)
covars$new_IMD_group <- as.factor(covars$new_IMD_group)
covars$comorbidity_category <- as.factor(covars$comorbidity_category)

combined_data <- cbind(covars, proteins)

# 4. Data Splitting for Model Training and Testing
set.seed(123)
split_index <- createDataPartition(combined_data$type, p = 0.7, list = FALSE)
train_data <- combined_data[split_index, ]
test_data <- combined_data[-split_index, ]

# 5. Data Scaling
## Scaling Selected Columns in Training and Testing Datasets
common_cols <- intersect(colnames(train_data), colnames(proteins))
additional_cols <- c('Age_at_recruitment', "BMI")
selected_cols <- c(common_cols, additional_cols)

train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)

# 6. Model Preparation and Cross-Validation
## Preparing Data for LASSO
train_proteins <- train_data[, intersect(colnames(train_data), colnames(proteins)), drop = FALSE]
train_covars <- train_data[, intersect(colnames(train_data), colnames(covars)), drop = FALSE]
y <- ifelse(train_covars$type == "CCB only", 1, 0)

confounders <- c('Age_at_recruitment', "Sex", "comorbidity_category", "Smoking_status", "BMI")
x <- cbind(train_proteins, train_covars[, confounders])
x <- model.matrix(~ . - 1, data = x)

penalty_factors <- rep(1, ncol(x))
penalty_factors[(ncol(x) - 8):ncol(x)] <- 0

## Cross-Validation for Model Selection
set.seed(1234)
t0 <- Sys.time()
mymodel <- cv.glmnet(x = x, y = y, penalty.factor = penalty_factors, type.measure = "auc", family = "binomial")
t1 <- Sys.time()
print(t1 - t0)

## Plotting Cross-Validation Results
plot(mymodel)

## Selected Variables Based on LASSO
beta_lasso <- coef(mymodel, s = "lambda.1se")[2:(ncol(proteins) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

# 7. Stability Selection and Feature Selection
## Running Stability Selection
t0 <- Sys.time()
out <- VariableSelection(xdata = x, ydata = y, verbose = FALSE, penalty.factor = penalty_factors, family = "binomial", n_cat = 3, k = 1000, pi_list = seq(0.5, 0.9, by = 0.01))
t1 <- Sys.time()
print(t1 - t0)

## Plotting Selection Proportions
CalibrationPlot(out)
hat_params <- Argmax(out)
selprop <- SelectionProportions(out)

par(mar = c(10, 5, 1, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "", ylab = "Selection Proportion", xaxt = "n", col = ifelse(selprop >= hat_params[2], "blue", "grey80"), cex.lab = 1)
abline(h = hat_params[2], lty = 2, col = "red")

threshold_indices <- which(selprop >= hat_params[2])
axis_labels <- names(selprop)[threshold_indices]

for (i in seq_along(threshold_indices)) {
  axis(side = 1, at = threshold_indices[i], labels = axis_labels[i], las = 2, col = "blue", col.axis = "blue", cex.axis = 0.7)
}

selected_protein <- names(selprop)[selprop >= 0.62]
selected_protein <- as.data.frame(selected_protein)

# 8. Logistic Regression Model Training and Evaluation
## Model Fitting and AUC Calculation on Training Data
train_data$type <- as.factor(train_data$type)
test_data$type <- as.factor(test_data$type)
train_data$type = as.factor(ifelse(as.character(train_data$type) == "CCB only", yes = 1, no = 0))
train_data$type <- relevel(train_data$type, ref = "0")
test_data$type = as.factor(ifelse(as.character(test_data$type) == "CCB only", yes = 1, no = 0))
test_data$type <- relevel(test_data$type, ref = "0")

lr.fit <- glm(type ~ CEACAM21 + CKMT1A_CKMT1B + CLUL1 + ECE1 + IL1R1 + IL7R + LEFTY2 + PDGFRB + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = train_data, family = "binomial")
summary(lr.fit)

## Model Prediction and ROC Curve for Training Data
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

## Model Fitting and AUC Calculation on Testing Data
lr.fit.test <- glm(type ~ CEACAM21 + CKMT1A_CKMT1B + CLUL1 + ECE1 + IL1R1 + IL7R + LEFTY2 + PDGFRB + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = test_data, family = "binomial")
summary(lr.fit.test)

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

# 9. Combined AUC Plot for Training and Testing Data
plot(train.perf.obj, col = "red", main = "", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.07, labels = paste0("Train AUC = ", round(as.numeric(train.perf.auc@y.values), 3)), col = "red", cex = 0.7)
plot(test.perf.obj, col = "blue", add = TRUE)
text(0.9, 0.0, labels = paste0("Test AUC = ", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

# 10. Repeated Training and Testing for Average AUC Calculation
auc_lr_total <- NULL
for (i in 1:100) {
  set.seed(6581 + i)
  split_index <- createDataPartition(combined_data$type, p = 0.7, list = FALSE)
  train_data <- combined_data[split_index, ]
  test_data <- combined_data[-split_index, ]
  
  train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
  test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)
  
  test_data$type = as.factor(ifelse(as.character(test_data$type) == "CCB only", yes = 1, no = 0))
  test_data$type <- relevel(test_data$type, ref = "0")
  
  lr.fit.test <- glm(type ~ CEACAM21 + CKMT1A_CKMT1B + CLUL1 + ECE1 + IL1R1 + IL7R + LEFTY2 + PDGFRB + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = test_data, family = "binomial")
  
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

# 11. Combined Training and Testing AUC Plot with Repeated Splits
auc_lr_train_total <- NULL
auc_lr_test_total <- NULL

for (i in 1:100) {
  set.seed(6581 + i)
  split_index <- createDataPartition(combined_data$type, p = 0.7, list = FALSE)
  train_data <- combined_data[split_index, ]
  test_data <- combined_data[-split_index, ]
  
  train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
  test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)
  
  train_data$type = as.factor(ifelse(as.character(train_data$type) == "CCB only", yes = 1, no = 0))
  test_data$type = as.factor(ifelse(as.character(test_data$type) == "CCB only", yes = 1, no = 0))
  
  lr.fit <- glm(type ~ CEACAM21 + CKMT1A_CKMT1B + CLUL1 + ECE1 + IL1R1 + IL7R + LEFTY2 + PDGFRB + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = train_data, family = "binomial")
  summary(lr.fit)
  
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
  
  lr.fit.test <- glm(type ~ CEACAM21 + CKMT1A_CKMT1B + CLUL1 + ECE1 + IL1R1 + IL7R + LEFTY2 + PDGFRB + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = test_data, family = "binomial")
  
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

# 12. Odds Ratio Calculation and Forest Plot Generation
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

forest_plot

# 13. Reordering and Renaming Variables for Forest Plot
reorder_rename_chr <- function(data, new_order, new_names) {
  for (i in seq_along(data)) {
    if (is.character(data[[i]]$Variable)) {
      data[[i]]$Variable <- factor(data[[i]]$Variable, levels = new_order, labels = new_names)
    }
  }
  return(data)
}

new_order <- c("Age_at_recruitment", "SexMale", "BMI", "Smoking_statusNever", "Smoking_statusPrevious", "comorbidity_category1", "comorbidity_category2", "comorbidity_category3+", "CEACAM21", "CKMT1A_CKMT1B", "CLUL1", "ECE1", "IL1R1", "IL7R", "LEFTY2", "PDGFRB")
new_names <- c("Age", "Male [ref = Female]", "BMI", "Smoking Status = Never [ref = Current]", "Smoking Status = Previous [ref = Current]", "Comorbidity Category 1 [ref = 0]", "Comorbidity Category 2 [ref = 0]", "Comorbidity Category +3 [ref = 0]", "CEACAM21", "CKMT1A_CKMT1B", "CLUL1", "ECE1", "IL1R1", "IL7R", "LEFTY2", "PDGFRB")

plot_data_list <- list(plot_data1, plot_data2)
modified_data <- reorder_rename_chr(plot_data_list, new_order, new_names)

plot_data1 <- modified_data[[1]]
plot_data2 <- modified_data[[2]]

# 14. Forest Plot Function and Plot Generation
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

remove_rows_by_rowname <- function(data, rowname_to_remove) {
  row_names <- rownames(data)
  data <- data[!(row_names %in% rowname_to_remove), , drop = FALSE]
  rownames(data) <- NULL
  return(data)
}

rowname_to_remove <- "(Intercept)"
plot_data1 <- remove_rows_by_rowname(plot_data1, rowname_to_remove)
plot_data2 <- remove_rows_by_rowname(plot_data2, rowname_to_remove)

plot_model1 <- create_forest_plot(plot_data1, "Training")
plot_model2 <- create_forest_plot(plot_data2, "Test")

combined_plot <- plot_model1 + labs(x = "Odds Ratio") | plot_model2 + labs(x = "Odds Ratio")
print(combined_plot)

# 15. Protein Expression Analysis: Upregulated and Downregulated Proteins
covars4 <- covars
protein_names <- selected_protein$selected_protein
MV_Sign_CCB <- proteins[, colnames(proteins) %in% protein_names]

MV_Sign_CCB <- rownames_to_column(MV_Sign_CCB, var = "EID")
covars4 <- rownames_to_column(covars4, var = "EID")
MV_Sign_CCB <- merge(MV_Sign_CCB, covars4[, c("EID", "type3")], by = "EID")
names(MV_Sign_CCB)[names(MV_Sign_CCB) == "type3"] <- "Treatment_Status"

summary(as.factor(MV_Sign_CCB$Treatment_Status))

protein_data_long <- melt(MV_Sign_CCB, id.vars = c("EID", "Treatment_Status"))
mean_values <- aggregate(value ~ variable + Treatment_Status, protein_data_long, mean, na.rm = TRUE)

step1_ccb_mean <- mean_values[mean_values$Treatment_Status == "CCB only", ]
untreated_mean <- mean_values[mean_values$Treatment_Status == "Untreated_CCB", ]
mean_diff <- merge(step1_ccb_mean, untreated_mean, by = "variable", suffixes = c("_Step1_CCB", "_Untreated"), all = TRUE)
mean_diff$Mean_Difference <- mean_diff$value_Step1_CCB - mean_diff$value_Untreated
mean_diff <- mean_diff[order(mean_diff$Mean_Difference, decreasing = TRUE), ]

top_upregulated <- head(mean_diff, 10)
top_downregulated <- tail(mean_diff, 10)

print("Top Upregulated Proteins in Step1 CCB:")
print(top_upregulated)

print("Top Downregulated Proteins in Step1 CCB:")
print(top_downregulated)

threshold <- 0.0
mean_diff$Regulation <- ifelse(mean_diff$Mean_Difference > threshold, "Upregulated", ifelse(mean_diff$Mean_Difference < -threshold, "Downregulated", "Unchanged"))

MV_CCB_Foldchange_Results <- as.data.frame(mean_diff)
write.table(MV_CCB_Foldchange_Results, "MV_CCB_Foldchange_Results.csv", sep = ",", row.names = FALSE)

# 16. Heatmap Visualization of Protein Expression
protein_data_long$value <- as.numeric(protein_data_long$value)
ggplot(mean_values, aes(x = Treatment_Status, y = variable, fill = value)) +
  geom_tile() +
  labs(x = "Treatment Status", y = "", fill = "NPX Value") +
  scale_fill_viridis_c() +
  scale_x_discrete(labels = c("CCB only" = "Step1 CCBs", "Untreated_CCB" = "Untreated-CCBs")) +
  theme_minimal()
