# Set working directory
setwd("/rds/general/user/sas123/home/TDS_Project/Results_Final_Project/260424_AIM2")
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
suppressPackageStartupMessages(library(patchwork))

# 2. Data Loading and Cleaning
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

# 3. Data Subsetting
subset_CCB_R <- covars$Response_Type %in% c("CCB_sufficient", "CCB_insufficient")
covars <- covars[subset_CCB_R, ]
row_names_subset <- rownames(covars)
proteins <- proteins[row_names_subset, ]

# Combine covariates and proteins data
combined_data <- cbind(covars, proteins)

# 4. Data Splitting
set.seed(123)
split_index <- createDataPartition(combined_data$Response_Type, p = 0.7, list = FALSE)
train_data <- combined_data[split_index, ]
test_data <- combined_data[-split_index, ]

# 5. Data Scaling
common_cols <- intersect(colnames(train_data), colnames(proteins))
additional_cols <- c("Age_at_recruitment", "BMI")
selected_cols <- c(common_cols, additional_cols)

train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)

# 6. Data Preparation for LASSO and Stability Selection
y <- ifelse(train_covars$Response_Type == "CCB_insufficient", 1, 0)
confounders <- c('Age_at_recruitment', "Sex", "comorbidity_category", "Smoking_status", "BMI")

x <- cbind(train_proteins, train_covars[, confounders])
x <- model.matrix(~ . - 1, data = x)

penalty_factors <- rep(1, ncol(x))
penalty_factors[(ncol(x) - 8):ncol(x)] <- 0
length(penalty_factors) == ncol(x)

# 7. LASSO Cross-validation
set.seed(1234)
t0 <- Sys.time()
mymodel <- cv.glmnet(x = x, y = y, penalty.factor = penalty_factors, type.measure = "deviance", family = "binomial", nfolds = 5)
t1 <- Sys.time()
print(t1 - t0)
plot(mymodel)

# Extracting selected variables
beta_lasso <- coef(mymodel, s = "lambda.min")[2:(ncol(proteins) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

# 8. Running Stability Selection
t0 <- Sys.time()
out <- VariableSelection(xdata = x, ydata = y, verbose = FALSE, penalty.factor = penalty_factors, family = "binomial", n_cat = 2, k = 1000, pi_list = seq(0.2, 0.6, by = 0.01))
t1 <- Sys.time()
print(t1 - t0)
CalibrationPlot(out)
hat_params <- Argmax(out)
print(hat_params)
selprop <- SelectionProportions(out)

# Visualization of selection proportions
par(mar = c(10, 5, 1, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "", ylab = "Selection Proportion", xaxt = "n", col = ifelse(selprop >= hat_params[2], "blue", "grey80"), cex.lab = 1)
abline(h = hat_params[2], lty = 2, col = "red")

threshold_indices <- which(selprop >= hat_params[2])
axis_labels <- names(selprop)[threshold_indices]
for (i in seq_along(threshold_indices)) {
  axis(side = 1, at = threshold_indices[i], labels = axis_labels[i], las = 2, col = "blue", col.axis = "blue", cex.axis = 0.7)
}

selected_protein <- names(selprop)[selprop >= 0.2]
print(selected_protein)
selected_protein <- as.data.frame(selected_protein)

# 9. Fold Change Calculation
covars_dem <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/covars_2.csv")
proteins_dem <- read.csv("/rds/general/user/yw423/projects/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/protein_data.csv")

covars_dem <- covars_dem[, -1, drop = FALSE]
proteins_dem <- proteins_dem[, -1, drop = FALSE]

subset_CCB_R <- covars_dem$Response_Type %in% c("CCB_sufficient", "CCB_insufficient")
covars_dem <- covars_dem[subset_CCB_R, ]
row_names_subset <- rownames(covars_dem)
proteins_dem <- proteins_dem[row_names_subset, ]

# Subset column names from proteins_dem based on selected_protein
protein_names <- selected_protein$selected_protein
MV_Sign_CCB <- proteins_dem[, c("EID", protein_names), drop = FALSE]
MV_Sign_CCB <- merge(MV_Sign_CCB, covars_dem[, c("EID", "Response_Type")], by = "EID")

# Reshape data to long format
protein_data_long <- melt(MV_Sign_CCB, id.vars = c("EID", "Response_Type"))

# Calculate the mean NPX value for each protein in each treatment group
mean_values <- aggregate(value ~ variable + Response_Type, protein_data_long, mean, na.rm = TRUE)

# Calculate the difference in mean NPX values
CCB_insufficient_mean <- mean_values[mean_values$Response_Type == "CCB_insufficient", ]
CCB_sufficient_mean <- mean_values[mean_values$Response_Type == "CCB_sufficient", ]
mean_diff <- merge(CCB_insufficient_mean, CCB_sufficient_mean, by = "variable", suffixes = c("_CCB_insufficient", "_CCB_sufficient"), all = TRUE)
mean_diff$Mean_Difference <- mean_diff$value_CCB_insufficient - mean_diff$value_CCB_sufficient

# Sort proteins and print top upregulated and downregulated proteins
mean_diff <- mean_diff[order(mean_diff$Mean_Difference, decreasing = TRUE), ]
top_upregulated <- head(mean_diff, 10)
top_downregulated <- tail(mean_diff, 10)
print("Top Upregulated Proteins in CCB_insufficient:")
print(top_upregulated)
print("Top Downregulated Proteins in CCB_insufficient:")
print(top_downregulated)

# Add regulation status
threshold <- 0.0
mean_diff$Regulation <- ifelse(mean_diff$Mean_Difference > threshold, "Upregulated", ifelse(mean_diff$Mean_Difference < -threshold, "Downregulated", "Unchanged"))

# Save results
MV_CCB_Foldchange_Results <- as.data.frame(mean_diff)
write.table(MV_CCB_Foldchange_Results, "MV_CCB_Foldchange_Results_Aim2.csv", sep = ",", row.names = FALSE)

# 10. Heatmap for Visualization
ggplot(mean_values, aes(x = Response_Type, y = variable, fill = value)) +
  geom_tile() +
  labs(x = "Medication Sufficiency", y = "Protein", fill = "NPX Value") +
  scale_fill_viridis_c() +
  scale_x_discrete(labels = c("CCB_insufficient" = "CCBs Insufficient", "CCB_sufficient" = "CCBs Sufficient")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# 11. Logistic Regression for Prediction

# Variable types
train_data$Response_Type <- as.factor(train_data$Response_Type)
test_data$Response_Type <- as.factor(test_data$Response_Type)

# Set the reference level for Response_Type
train_data$Response_Type <- relevel(train_data$Response_Type, ref = "CCB_sufficient")
test_data$Response_Type <- relevel(test_data$Response_Type, ref = "CCB_sufficient")
train_data$Response_Type <- as.factor(ifelse(as.character(train_data$Response_Type) == "CCB_insufficient", yes = 1, no = 0))
test_data$Response_Type <- as.factor(ifelse(as.character(test_data$Response_Type) == "CCB_insufficient", yes = 1, no = 0))

## Train data ##
lr.fit <- glm(Response_Type ~ PNLIPRP2 + CCL28 + DDX58 + NTF3 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, 
              data = train_data, family = "binomial")
summary(lr.fit)

# Predict probabilities and plot ROC for training data
lr.train_pred <- predict(lr.fit, type = "response", newdata = train_data)
train_data <- cbind(train_data, probs = lr.train_pred)

train.pred.obj <- prediction(train_data$probs, train_data$Response_Type)
train.perf.obj <- performance(train.pred.obj, "tpr", "fpr")
train.perf.auc <- performance(train.pred.obj, "auc")
plot(train.perf.obj, col = "blue", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.05, labels = paste0("AUC=", round(as.numeric(train.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

# Confusion matrix for training data
train_predicted_class <- factor(ifelse(lr.train_pred > 0.5, "1", "0"), levels = levels(train_data$Response_Type))
train_conf_matrix <- confusionMatrix(train_predicted_class, train_data$Response_Type)
print(train_conf_matrix)

## Test Data ##
lr.fit.test <- glm(Response_Type ~ PNLIPRP2 + CCL28 + DDX58 + NTF3 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, 
                   data = test_data, family = "binomial")
summary(lr.fit.test)

# Predict probabilities and plot ROC for testing data
lr.test_pred <- predict(lr.fit.test, type = "response", newdata = test_data)
test_data <- cbind(test_data, probs = lr.test_pred)

test.pred.obj <- prediction(test_data$probs, test_data$Response_Type)
test.perf.obj <- performance(test.pred.obj, "tpr", "fpr")
test.perf.auc <- performance(test.pred.obj, "auc")
plot(test.perf.obj, col = "blue", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.05, labels = paste0("AUC=", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

# Confusion matrix for testing data
test_predicted_class <- factor(ifelse(lr.test_pred > 0.5, "1", "0"), levels = levels(test_data$Response_Type))
test_conf_matrix <- confusionMatrix(test_predicted_class, test_data$Response_Type)
print(test_conf_matrix)

# 12. Combined AUC Plot
plot(train.perf.obj, col = "red", main = "", xlab = "1 - Specificity", ylab = "Sensitivity")
abline(a = 0, b = 1, lty = "dashed", col = "gray")
text(0.9, 0.07, labels = paste0("Train AUC = ", round(as.numeric(train.perf.auc@y.values), 3)), col = "red", cex = 0.7)
plot(test.perf.obj, col = "blue", add = TRUE)
text(0.9, 0.0, labels = paste0("Test AUC = ", round(as.numeric(test.perf.auc@y.values), 3)), col = "blue", cex = 0.7)

# 13. Repeated AUC Calculation with Multiple Seeds
auc_lr_total <- NULL
for (i in 1:100) {
  set.seed(6581 + i)
  split_index <- createDataPartition(combined_data$Response_Type, p = 0.7, list = FALSE)
  train_data <- combined_data[split_index, ]
  test_data <- combined_data[-split_index, ]
  train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
  test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)
  
  train_data$Response_Type <- as.factor(ifelse(as.character(train_data$Response_Type) == "CCB_insufficient", yes = 1, no = 0))
  train_data$Response_Type <- relevel(train_data$Response_Type, ref = "0")
  
  test_data$Response_Type <- as.factor(ifelse(as.character(test_data$Response_Type) == "CCB_insufficient", yes = 1, no = 0))
  test_data$Response_Type <- relevel(test_data$Response_Type, ref = "0")
  
  lr.fit.test <- glm(Response_Type ~ PNLIPRP2 + CCL28 + DDX58 + NTF3 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, 
                     data = test_data, family = "binomial")
  
  lr.test_pred <- predict(lr.fit.test, type = "response", newdata = test_data)
  test_data <- cbind(test_data, probs = lr.test_pred)
  
  test.pred.obj <- prediction(test_data$probs, test_data$Response_Type)
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

# 14. Combined AUC Plot for Train and Test
auc_lr_train_total <- NULL
auc_lr_test_total <- NULL
for (i in 1:100) {
  set.seed(6581 + i)
  split_index <- createDataPartition(combined_data$Response_Type, p = 0.7, list = FALSE)
  train_data <- combined_data[split_index, ]
  test_data <- combined_data[-split_index, ]
  train_data[, selected_cols] <- apply(train_data[, selected_cols], 2, scale)
  test_data[, selected_cols] <- apply(test_data[, selected_cols], 2, scale)
  
  train_data$Response_Type <- as.factor(ifelse(as.character(train_data$Response_Type) == "CCB_insufficient", yes = 1, no = 0))
  test_data$Response_Type <- as.factor(ifelse(as.character(test_data$Response_Type) == "CCB_insufficient", yes = 1, no = 0))
  
  lr.fit_train <- glm(Response_Type ~ PNLIPRP2 + CCL28 + DDX58 + NTF3 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, 
                      data = train_data, family = "binomial")
  
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
  
  lr.fit_test <- glm(Response_Type ~ PNLIPRP2 + CCL28 + DDX58 + NTF3 + Age_at_recruitment + Sex + BMI + Smoking_status + comorbidity_category, 
                     data = test_data, family = "binomial")
  
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

# 15. Forest Plot for Coefficients
extract_coefficients <- function(model_summary) {
  coefficients <- model_summary$coefficients[, 1]
  std_errors <- model_summary$coefficients[, 2]
  odds_ratios <- exp(coefficients)
  ci_lower <- exp(coefficients - 1.96 * std_errors)
  ci_upper <- exp(coefficients + 1.96 * std_errors)
  variables <- rownames(model_summary$coefficients)
  
  df <- data.frame(Variable = variables, Coefficient = coefficients, OR = odds_ratios, CI_lower = ci_lower, CI_upper = ci_upper)
  return(df)
}

plot_data1 <- extract_coefficients(summary(lr.fit))
plot_data1$Model <- "Training"
plot_data2 <- extract_coefficients(summary(lr.fit.test))
plot_data2$Model <- "Test"

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
  theme(axis.text.y = element_text(size = 8)) +
  xlim(0.1, 10)

print(forest_plot)

# 16. Reordering and Renaming for Forest Plot
reorder_rename_chr <- function(data, new_order, new_names) {
  for (i in seq_along(data)) {
    if (is.character(data[[i]]$Variable)) {
      data[[i]]$Variable <- factor(data[[i]]$Variable, levels = new_order, labels = new_names)
    }
  }
  return(data)
}

new_order <- c("Age_at_recruitment", "SexMale", "BMI", "Smoking_statusNever", "Smoking_statusPrevious", "comorbidity_category1", "comorbidity_category2", "comorbidity_category3+", "CCL28", "DDX58", "NTF3", "PNLIPRP2")
new_names <- c("Age", "Male [ref = Female]", "BMI", "Smoking Status = Never [ref = Current]", "Smoking Status = Previous [ref = Current]", "Comorbidity Category 1 [ref = 0]", "Comorbidity Category 2 [ref = 0]", "Comorbidity Category +3 [ref = 0]", "CCL28", "DDX58", "NTF3", "PNLIPRP2")

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

# 17. Matching and Stability Selection
m.out <- matchit(Response_Type ~ Age_at_recruitment + BMI, data = train_data, method = "nearest", exact = ~Sex, distance = "euclidean", ratio = 3)
matched_samples_train <- match.data(m.out)
matched_samples_train$Response_Type <- as.factor(matched_samples_train$Response_Type)

matched_samples_train[, selected_cols] <- apply(matched_samples_train[, selected_cols], 2, scale)
train_proteins <- matched_samples_train[, intersect(colnames(matched_samples_train), colnames(proteins)), drop = FALSE]
train_covars <- matched_samples_train[, intersect(colnames(matched_samples_train), colnames(covars)), drop = FALSE]

y <- train_covars$Response_Type
x <- cbind(train_proteins, train_covars[, confounders])
x <- model.matrix(~ . - 1, data = x)

penalty_factors <- rep(1, ncol(x))
penalty_factors[(ncol(x) - 8):ncol(x)] <- 0

set.seed(1234)
t0 <- Sys.time()
mymodel <- cv.glmnet(x = x, y = y, penalty.factor = penalty_factors, type.measure = "deviance", nfolds = 5, family = "binomial")
t1 <- Sys.time()
print(t1 - t0)
plot(mymodel)

beta_lasso <- coef(mymodel, s = "lambda.min")[2:(ncol(proteins) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

# Stability selection
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
