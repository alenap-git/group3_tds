setwd("/rds/general/user/sas123/home/TDS_Project/Results_Final_Project/260424_AIM2")
getwd()

# 1. Loading Libraries and Cleaning the Workspace
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
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(patchwork))

# 2. Loading and Preprocessing Data
covars <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/covars_2.csv")
proteins <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/protein_data.csv")

covars <- covars[, -1, drop = FALSE]
rownames(covars) <- covars[, 1]
covars <- covars[, -1, drop = FALSE]

proteins <- proteins[, -1, drop = FALSE]
rownames(proteins) <- proteins[, 1]
proteins <- proteins[, -1, drop = FALSE]

covars$Sex <- as.factor(covars$Sex)
covars$modified_ethnicity <- as.factor(covars$modified_ethnicity)
covars$Alcohol_intake_recode <- as.factor(covars$Alcohol_intake_recode)
covars$Smoking_status <- as.factor(covars$Smoking_status)
covars$new_IMD_group <- as.factor(covars$new_IMD_group)
covars$comorbidity_category <- as.factor(covars$comorbidity_category)

summary(as.factor(covars$Response_Type))

subset_ACE_R <- covars$Response_Type %in% c("ACE/ARBs_sufficient", "ACE/ARBs_insufficient")
covars <- covars[subset_ACE_R, ]

row_names_subset <- rownames(covars)
proteins <- proteins[row_names_subset, ]

combined_data <- cbind(covars, proteins)
summary(as.factor(combined_data$Response_Type))

# 3. Splitting Data into Training and Testing Sets
set.seed(123)
split_index <- createDataPartition(combined_data$Response_Type, p = 0.7, list = FALSE)
train_data <- combined_data[split_index, ]
test_data <- combined_data[-split_index, ]

summary(as.factor(train_data$Response_Type))
summary(as.factor(test_data$Response_Type))

# 4. Data Scaling and Preparation for Modeling
common_cols <- intersect(colnames(train_data), colnames(proteins))
additional_cols <- c('Age_at_recruitment', "BMI")
selected_cols <- c(common_cols, additional_cols)

train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)

train_proteins <- train_data[, intersect(colnames(train_data), colnames(proteins)), drop = FALSE]
train_covars <- train_data[, intersect(colnames(train_data), colnames(covars)), drop = FALSE]

y <- ifelse(train_covars$Response_Type == "ACE/ARBs_insufficient", 1, 0)
confounders <- c('Age_at_recruitment', "Sex", "comorbidity_category", "Smoking_status", "BMI")

contrasts(as.factor(train_covars$Response_Type))

x <- cbind(train_proteins, train_covars[, confounders])
x <- model.matrix(~ . - 1, data = x)

penalty_factors <- rep(1, ncol(x))
penalty_factors[(ncol(x) - 8):ncol(x)] <- 0

length(penalty_factors) == ncol(x)
v_penalty_factors <- tail(penalty_factors, 15)

# 5. Model Training and Cross-Validation
set.seed(1234)
t0 <- Sys.time()
mymodel <- cv.glmnet(x = x, y = y, penalty.factor = penalty_factors, type.measure = "auc", family = "binomial")
t1 <- Sys.time()
print(t1 - t0)

plot(mymodel)

# 6. Variable Selection and Stability Selection
beta_lasso <- coef(mymodel, s = "lambda.1se")[2:(ncol(proteins) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

t0 <- Sys.time()
out <- VariableSelection(xdata = x, ydata = y, verbose = FALSE, penalty.factor = penalty_factors, family = "binomial", n_cat = 2, k = 1000, pi_list = seq(0.4, 0.9, by = 0.01))
t1 <- Sys.time()
print(t1 - t0)

CalibrationPlot(out)
hat_params <- Argmax(out)
print(hat_params)

selprop <- SelectionProportions(out)

# 7. Visualization of Selection Proportions
par(mar = c(10, 5, 1, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "", ylab = "Selection Proportion", xaxt = "n", col = ifelse(selprop >= hat_params[2], "blue", "grey80"), cex.lab = 1)
abline(h = hat_params[2], lty = 2, col = "red")

threshold_indices <- which(selprop >= hat_params[2])
axis_labels <- names(selprop)[threshold_indices]

for (i in seq_along(threshold_indices)) {
  axis(side = 1, at = threshold_indices[i], labels = axis_labels[i], las = 2, col = "blue", col.axis = "blue", cex.axis = 0.7)
}

selected_protein <- names(selprop)[selprop >= 0.4]
print(selected_protein)
selected_protein <- as.data.frame(selected_protein)

# 8. Protein Expression Analysis
covars_dem <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/covars_2.csv")
proteins_dem <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/protein_data.csv")

covars_dem <- covars_dem[, -1, drop = FALSE]
proteins_dem <- proteins_dem[, -1, drop = FALSE]

subset_ACE_R <- covars_dem$Response_Type %in% c("ACE/ARBs_sufficient", "ACE/ARBs_insufficient")
covars_dem <- covars_dem[subset_ACE_R, ]

row_names_subset <- rownames(covars_dem)
proteins_dem <- proteins_dem[row_names_subset, ]

protein_names <- selected_protein$selected_protein
MV_Sign_ACE <- proteins_dem[, c("EID", protein_names), drop = FALSE]

MV_Sign_ACE <- merge(MV_Sign_ACE, covars_dem[, c("EID", "Response_Type")], by = "EID")

protein_data_long <- melt(MV_Sign_ACE, id.vars = c("EID", "Response_Type"))
mean_values <- aggregate(value ~ variable + Response_Type, protein_data_long, mean, na.rm = TRUE)

ACE_ARBs_insufficient_mean <- mean_values[mean_values$Response_Type == "ACE/ARBs_insufficient", ]
ACE_ARBs_sufficient_mean <- mean_values[mean_values$Response_Type == "ACE/ARBs_sufficient", ]

mean_diff <- merge(ACE_ARBs_insufficient_mean, ACE_ARBs_sufficient_mean, by = "variable", suffixes = c("_ACE_ARBs_insufficient", "_ACE_ARBs_sufficient"), all = TRUE)
mean_diff$Mean_Difference <- mean_diff$value_ACE_ARBs_insufficient - mean_diff$value_ACE_ARBs_sufficient
mean_diff <- mean_diff[order(mean_diff$Mean_Difference, decreasing = TRUE), ]

top_upregulated <- head(mean_diff, 10)
top_downregulated <- tail(mean_diff, 10)

print("Top Upregulated Proteins in Step1 ACE:")
print(top_upregulated)

print("Top Downregulated Proteins in Step1 ACE:")
print(top_downregulated)

threshold <- 0.0
mean_diff$Regulation <- ifelse(mean_diff$Mean_Difference > threshold, "Upregulated", ifelse(mean_diff$Mean_Difference < -threshold, "Downregulated", "Unchanged"))

head(mean_diff)
Aim2_MV_ACE_Foldchange_Results <- as.data.frame(mean_diff)
write.table(Aim2_MV_ACE_Foldchange_Results, "MV_ACE_Foldchange_Results_Aim2.csv", sep = ",", row.names = FALSE)

# 9. Heatmap Visualization of Protein Expression
ggplot(mean_values, aes(x = Response_Type, y = variable, fill = value)) +
  geom_tile() +
  labs(x = "Medication Sufficiency", y = "Protein", fill = "NPX Value") +
  scale_fill_viridis_c() +
  scale_x_discrete(labels = c("ACE/ARBs_insufficient" = "ACE/ARBs Insufficient", "ACE/ARBs_sufficient" = "ACE/ARBs Sufficient")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# 10. T-test Analysis for Protein Expression Differences
result_list <- lapply(unique(protein_data_long$variable), function(var) {
  df_subset <- protein_data_long[protein_data_long$variable == var, ]
  t_test_result <- t.test(value ~ Response_Type, data = df_subset)
  return(list(variable = var, t_test = t_test_result))
})

for (res in result_list) {
  cat("Variable:", res$variable, "\n")
  print(res$t_test)
  cat("\n")
}

# 11. Logistic Regression Modeling
train_data$Response_Type <- as.factor(train_data$Response_Type)
test_data$Response_Type <- as.factor(test_data$Response_Type)

train_data$Response_Type <- relevel(train_data$Response_Type, ref = "ACE/ARBs_sufficient")
test_data$Response_Type <- relevel(test_data$Response_Type, ref = "ACE/ARBs_sufficient")

train_data$Response_Type = as.factor(ifelse(as.character(train_data$Response_Type) == "ACE/ARBs_insufficient", yes = 1, no = 0))
test_data$Response_Type = as.factor(ifelse(as.character(test_data$Response_Type) == "ACE/ARBs_insufficient", yes = 1, no = 0))

lr.fit <- glm(Response_Type ~ TDGF1 + SERPINA9 + LTA + LILRB5 + IL1RAP + COL9A1 + CD200R1 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = train_data, family = "binomial")
summary(lr.fit)

# 12. ROC Curve for Training Data
lr.train_pred <- predict(lr.fit, type = "response", newdata = train_data)
train_data <- cbind(train_data, probs = lr.train_pred)

train.pred.obj <- prediction(train_data$probs, train_data$Response_Type)
train.perf.obj <- performance(train.pred.obj, "tpr", "fpr")
train.perf.auc <- performance(train.pred.obj, "auc")

plot(train.perf.obj, col = "blue", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.05, labels = paste0("AUC=", round(as.numeric(train.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

train_predicted_class <- factor(ifelse(lr.train_pred > 0.5, "1", "0"), levels = levels(train_data$Response_Type))
train_conf_matrix <- confusionMatrix(train_predicted_class, train_data$Response_Type)
print(train_conf_matrix)

# 13. ROC Curve for Testing Data
lr.fit.test <- glm(Response_Type ~ TDGF1 + SERPINA9 + LTA + LILRB5 + IL1RAP + COL9A1 + CD200R1 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = test_data, family = "binomial")
summary(lr.fit.test)

lr.test_pred <- predict(lr.fit.test, type = "response", newdata = test_data)
test_data <- cbind(test_data, probs = lr.test_pred)

test.pred.obj <- prediction(test_data$probs, test_data$Response_Type)
test.perf.obj <- performance(test.pred.obj, "tpr", "fpr")
test.perf.auc <- performance(test.pred.obj, "auc")

plot(test.perf.obj, col = "blue", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.05, labels = paste0("AUC=", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

test_predicted_class <- factor(ifelse(lr.test_pred > 0.5, "1", "0"), levels = levels(test_data$Response_Type))
test_conf_matrix <- confusionMatrix(test_predicted_class, test_data$Response_Type)
print(test_conf_matrix)

# 14. Combined AUC Plot for Training and Testing Data
plot(train.perf.obj, col = "red", main = "", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.07, labels = paste0("Train AUC = ", round(as.numeric(train.perf.auc@y.values), 3)), col = "red", cex = 0.7)
plot(test.perf.obj, col = "blue", add = TRUE)
text(0.9, 0.0, labels = paste0("Test AUC = ", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

# 15. Cross-Validation for Training and Testing Data
auc_lr_total <- NULL
for (i in 1:100) {
  set.seed(6581 + i)
  split_index <- createDataPartition(combined_data$Response_Type, p = 0.7, list = FALSE)
  train_data <- combined_data[split_index, ]
  test_data <- combined_data[-split_index, ]
  
  train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
  test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)
  
  train_data$Response_Type = as.factor(ifelse(as.character(train_data$Response_Type) == "ACE/ARBs_insufficient", yes = 1, no = 0))
  test_data$Response_Type = as.factor(ifelse(as.character(test_data$Response_Type) == "ACE/ARBs_insufficient", yes = 1, no = 0))
  
  lr.fit_train <- glm(Response_Type ~ TDGF1 + SERPINA9 + LTA + LILRB5 + IL1RAP + COL9A1 + CD200R1 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = train_data, family = "binomial")
  
  lr.train_pred <- predict(lr.fit_train, type = "response", newdata = train_data)
  train_data <- cbind(train_data, probs = lr.train_pred)
  
  train.pred.obj <- prediction(train_data$probs, train_data$Response_Type)
  train.perf.obj <- performance(train.pred.obj, "tpr", "fpr")
  train.perf.auc <- performance(train.pred.obj, "auc")
  
  auc_lr_train <- train.perf.auc@y.values[[1]]
  auc_lr_train_total <- c(auc_lr_train_total, auc_lr_train)
  
  if (i == 1) {
    plot(train.perf.obj, col = "red", xlab = "1 - Specificity", ylab = "Sensitivity")
  } else {
    plot(train.perf.obj, add = TRUE, col = "red")
  }
  
  lr.fit_test <- glm(Response_Type ~ TDGF1 + SERPINA9 + LTA + LILRB5 + IL1RAP + COL9A1 + CD200R1 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = test_data, family = "binomial")
  
  lr.test_pred <- predict(lr.fit_test, type = "response", newdata = test_data)
  test_data <- cbind(test_data, probs = lr.test_pred)
  
  test.pred.obj <- prediction(test_data$probs, test_data$Response_Type)
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

# 16. Forest Plot of Logistic Regression Coefficients
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
  facet_wrap(~ Model, scales = "free_x") +
  theme(axis.text.y = element_text(size = 8))

forest_plot <- forest_plot + xlim(0.1, 15)
print(forest_plot)

# 17. Custom Forest Plot with Variable Reordering
new_order <- c("Age_at_recruitment", "SexMale", "BMI", "Smoking_statusNever", "Smoking_statusPrevious", "comorbidity_category1", "comorbidity_category2", "comorbidity_category3+", "CD200R1", "COL9A1", "IL1RAP", "LILRB5", "LTA", "SERPINA9", "TDGF1")
new_names <- c("Age", "Male [ref = Female]", "BMI", "Smoking Status = Never [ref = Current]", "Smoking Status = Previous [ref = Current]", "Comorbidity Category 1 [ref = 0]", "Comorbidity Category 2 [ref = 0]", "Comorbidity Category +3 [ref = 0]", "CD200R1", "COL9A1", "IL1RAP", "LILRB5", "LTA", "SERPINA9", "TDGF1")

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

plot_data1 <- remove_rows_by_rowname(plot_data1, rowname_to_remove)
plot_data2 <- remove_rows_by_rowname(plot_data2, rowname_to_remove)

plot_model1 <- create_forest_plot(plot_data1, "Training")
plot_model2 <- create_forest_plot(plot_data2, "Test")

combined_plot <- plot_model1 + labs(x = "Odds Ratio") | plot_model2 + labs(x = "Odds Ratio")
print(combined_plot)

# 18. Matching and Logistic Regression
m.out <- matchit(Response_Type ~ Age_at_recruitment + BMI, data = train_data, method = "nearest", exact = ~Sex, distance = "euclidean", ratio = 3)
matched_samples_train <- match.data(m.out)
summary(as.factor(matched_samples_train$Response_Type))

matched_samples_train[, selected_cols] <- apply(matched_samples_train[, selected_cols], 2, scale)

train_proteins <- matched_samples_train[, intersect(colnames(matched_samples_train), colnames(proteins)), drop = FALSE]
train_covars <- matched_samples_train[, intersect(colnames(matched_samples_train), colnames(covars)), drop = FALSE]

y <- train_covars$Response_Type
confounders <- c('Age_at_recruitment', "Sex", "comorbidity_category", "Smoking_status", "BMI")
x <- cbind(train_proteins, train_covars[, confounders])
x <- model.matrix(~ . - 1, data = x)

penalty_factors <- rep(1, ncol(x))
penalty_factors[(ncol(x) - 8):ncol(x)] <- 0
length(penalty_factors) == ncol(x)
v_penalty_factors <- tail(penalty_factors, 15)

set.seed(1234)
t0 <- Sys.time()
mymodel <- cv.glmnet(x = x, y = y, penalty.factor = penalty_factors, type.measure = "auc", nfolds = 10, family = "binomial")
t1 <- Sys.time()
print(t1 - t0)

plot(mymodel)

beta_lasso <- coef(mymodel, s = "lambda.1se")[2:(ncol(proteins) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

t0 <- Sys.time()
out <- VariableSelection(xdata = x, ydata = y, verbose = FALSE, penalty.factor = penalty_factors, family = "binomial", n_cat = 3, pi_list = seq(0.5, 0.9, by = 0.01))
t1 <- Sys.time()
print(t1 - t0)

CalibrationPlot(out)
hat_params <- Argmax(out)
print(hat_params)

selprop <- SelectionProportions(out)

par(mar = c(10, 5, 1, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "", ylab = "Selection Proportion", xaxt = "n", col = ifelse(selprop >= hat_params[2], yes = "red", no = "grey"), cex.lab = 1.5)
abline(h = hat_params[2], lty = 2, col = "darkred")

for (i in 1:length(selprop)) {
  axis(side = 1, at = i, labels = names(selprop)[i], las = 2, col = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"), col.axis = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"), cex.axis = 0.5)
}

selected_protein <- names(selprop)[selprop >= 0.1]
print(selected_protein)

# 19. Logistic Regression with Matched Samples
matched_samples_train$Response_Type <- as.factor(matched_samples_train$Response_Type)
test_data$Response_Type <- as.factor(test_data$Response_Type)

summary(matched_samples_train$Response_Type)
matched_samples_train$Response_Type <- relevel(matched_samples_train$Response_Type, ref = "ACE/ARBs_sufficient")
test_data$Response_Type <- relevel(test_data$Response_Type, ref = "ACE/ARBs_sufficient")

lr.fit <- glm(Response_Type ~ LILRB5 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, data = matched_samples_train, family = "binomial")
summary(lr.fit)

lr.pred.probs <- predict(lr.fit, type = "response", newdata = matched_samples_train)
matched_samples_train <- cbind(matched_samples_train, probs = lr.pred.probs)

train.pred.obj <- prediction(matched_samples_train$probs, matched_samples_train$Response_Type)
train.perf.obj <- performance(train.pred.obj, "tpr", "fpr")
train.perf.auc <- performance(train.pred.obj, "auc")

plot(train.perf.obj, col = "blue")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.05, labels = paste0("AUC=", round(as.numeric(train.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

train_predicted_class <- factor(ifelse(lr.train_pred > 0.5, "ACE/ARBs_insufficient", "ACE/ARBs_sufficient"), levels = levels(train_data$Response_Type))
train_conf_matrix <- confusionMatrix(train_predicted_class, train_data$Response_Type)
print(train_conf_matrix)

logit_model_Pred <- predict(lr.fit, test_data, type = "response")
predicted_class <- factor(ifelse(logit_model_Pred > 0.5, "ACE/ARBs_insufficient", "ACE/ARBs_sufficient"), levels = levels(test_data$Response_Type))
test_conf_matrix <- confusionMatrix(test_predicted_class, test_data$Response_Type)
print(test_conf_matrix)

test.pred.obj <- prediction(logit_model_Pred, test_data$Response_Type)
test.perf.obj <- performance(test.pred.obj, "tpr", "fpr")
test.perf.auc <- performance(test.pred.obj, "auc")

plot(test.perf.obj, col = "blue")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.0, labels = paste0("AUC = ", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.8)

# 20. Combined AUC Plot for Matched Training and Testing Data
plot(train.perf.obj, col = "red", main = "ROC Curve - Training (matched) vs Test Data")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.05, labels = paste0("Train AUC = ", round(as.numeric(train.perf.auc@y.values), 3)), col = "red", cex = 0.7)
plot(test.perf.obj, col = "blue", add = TRUE)
text(0.9, 0.0, labels = paste0("Test AUC = ", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

# 21. Additional Model Evaluation Metrics
logit_model_Pred <- predict(lr.fit, test_data, type = "response")
predicted_class <- ifelse(logit_model_Pred > 0.5, 1, 0)
pred.obj = prediction(logit_model_Pred, test_data$Response_Type)
perf.obj = performance(pred.obj, "tpr", "fpr")
perf.auc <- performance(pred.obj, "auc")

plot(perf.obj, col = "blue")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.0, labels = paste0("AUC=", round(as.numeric(perf.auc@y.values), 3)), col = "blue", cex = 0.8)

conf_matrix <- table(predicted_class, test_data$Response_Type)
conf_matrix

accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
specificity <- conf_matrix[1, 1] / sum(conf_matrix[1, ])
precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
recall <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
f1_score <- 2 * (precision * recall) / (precision + recall)

print(conf_matrix)
cat("\nAccuracy:", accuracy)
cat("\nSensitivity:", sensitivity)
cat("\nSpecificity:", specificity)
cat("\nPrecision:", precision)
cat("\nF1 Score:", f1_score)
