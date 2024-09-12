

### Setting the working directory
wd = "/rds/general/project/hda_23-24/live/TDS/group3/TDS_group3/Dataset"  
setwd(wd)  
getwd()  

### Loading and exploring the data
Data_1 = read.csv("clean_data.csv")  
dim(Data_1)  
head(Data_1) 
str(Data_1)  
print(rownames(Data_1)[1:10])  
levels(as.factor(Data_1$Comorbidities))  # Display the unique levels of the 'Comorbidities' variable

### checking the missing values
any(is.na(Data_1)) 
sum(is.na(Data_1))  
colSums(is.na(Data_1))  


### Cleaning and recoding 'Baseline_medication' as a factor
Data_1$Baseline_medication <- factor(Data_1$Baseline_medication)  
Data_1_indiv = Data_1[!duplicated(Data_1$EID),]  # Remove duplicate entries based on 'EID'


### Categorizing Diseases
categorize_disease <- function(disease_name) {
  if (disease_name %in% c("Arteries, arterioles, and capillary diseases", "Cardiac arrhythmias", "Cerebrovascular disease",
                          "Chronic rheumatic heart diseases", "Coagulopathy", "Hypertensive diseases", 
                          "Ischaemic heart diseases", "Other circulatory system disorders", "Pulmonary circulation disorders",
                          "Other heart failure", "Veins, lymphatic vessels and lymph node diseases")) {
    return("Cardiovascular and Vascular Disease")
  } else if (disease_name %in% c("Diabetes mellitus", "Disorders of other endocrine glands", "Other metabolic disorders",
                                 "Nutritional deficiencies", "Obesity", "Thyroid gland disorders")) {
    return("Metabolic Disorders")
  } else {
    return("Others")
  }
}

Data_1_indiv$Category <- sapply(Data_1_indiv$Disease, categorize_disease)  # Apply the categorization function to the 'Disease' column

### Summarizing by 'Baseline_medication' and 'Category'
library(dplyr)
result <- Data_1_indiv %>%
  group_by(Baseline_medication, Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Proportion = paste0(round((Count / sum(Count)) * 100, 2), "%"))  # Calculate the count and proportion for each group

print(result)  # Print the summary table

### Recoding 'Comorbidities' Variable
comorbidities_recode_rules_1 <- c(
  "Arteries, arterioles, and capillary diseases" = "Cardiovascular and vascular disease",
  "Cardiac arrhythmias" = "Cardiovascular and vascular disease",
  "Cerebrovascular disease" = "Cardiovascular and vascular disease",
  "Chronic rheumatic heart diseases" = "Cardiovascular and vascular disease",
  "Coagulopathy" = "Cardiovascular and vascular disease",
  "Hypertensive diseases" = "Cardiovascular and vascular disease",
  "Ischaemic heart diseases" = "Cardiovascular and vascular disease",
  "Other circulatory system disorders" = "Cardiovascular and vascular disease",
  "Pulmonary circulation disorders" = "Cardiovascular and vascular disease",
  "Other heart failure" = "Cardiovascular and vascular disease",
  "Veins, lymphatic vessels and lymph node diseases" = "Cardiovascular and vascular disease",
  "Kidney Disease" = "Renal and Urinary System Disorders",
  "Other urinary system diseases" = "Renal and Urinary System Disorders",
  "Diabetes mellitus" = "Metabolic disorders",
  "Disorders of other endocrine glands" = "Metabolic disorders",
  "Other metabolic disorders" = "Metabolic disorders",
  "Nutritional deficiencies" = "Metabolic disorders",
  "Obesity" = "Metabolic disorders",
  "Thyroid gland disorders" = "Metabolic disorders",
  "Blood loss anemia" = "Hematologic Conditions",
  "Deficiency anemia" = "Hematologic Conditions",
  "Other anemia" = "Hematologic Conditions",
  "Other diseases of blood and blood-forming organs" = "Hematologic Conditions",
  "HIV/AIDS" = "Immune System and Inflammatory Diseases",
  "Immune mechanism disorders" = "Immune System and Inflammatory Diseases",
  "Inflammatory diseases of CNS" = "Immune System and Inflammatory Diseases",
  "Rheumatoid arthritis" = "Immune System and Inflammatory Diseases",
  "Other arthritis" = "Immune System and Inflammatory Diseases",
  "Metastatic cancer" = "Cancer and Neoplastic Diseases",
  "Other neoplasms" = "Cancer and Neoplastic Diseases",
  "Depression" = "Neurological and Mental Health Disorders",
  "Mental disorder due to alcohol" = "Neurological and Mental Health Disorders",
  "Nerve, nerve root and plexus disorders" = "Neurological and Mental Health Disorders",
  "Other nervous system disorders" = "Neurological and Mental Health Disorders",
  "Other mental disorders" = "Neurological and Mental Health Disorders",
  "Fluid and electrolyte disorders" = "Others",
  "Other digestive system diseases" = "Others",
  "Other infectious and parasitic diseases" = "Others",
  "Other musculoskeletal diseases" = "Others",
  "Peptic ulcer disease" = "Others",
  "Respiratory system diseases" = "Others",
  "Other" = "Others"
)

# Recode 'Comorbidities' based on the defined rules

Data_1_indiv$Comorbidities_recode <- recode(Data_1_indiv$Comorbidities, !!!comorbidities_recode_rules_1)  

### Checking: 'Baseline_medication' and Recoded 'Comorbidities'
analysis_results <- lapply(split(Data_1_indiv, Data_1_indiv$Baseline_medication), function(subset) {
  tmp <- table(subset$Comorbidities_recode)
  data.frame(Comorbidity = names(tmp), Count = as.integer(tmp), Percentage = round((as.integer(tmp) / sum(tmp) * 100), 2))
})

print(analysis_results)  # Print the analysis results

### Recoding 'Comorbidities' Categories at boarder range 
Data_1_indiv <- Data_1_indiv %>%
  mutate(Comorbidities_cato = case_when(
    Full_Comorbidities %in% c("Ischaemic heart diseases", "Cerebrovascular disease", 
                              "Veins, lymphatic vessels and lymph node diseases", 
                              "Chronic rheumatic heart diseases", "Pulmonary circulation disorders", 
                              "Cardiovascular diseases", "Other circulatory system disorders", 
                              "Diabetes mellitus", "Cardiac arrhythmias", "Coagulopathy", 
                              "Arteries, arterioles, and capillary diseases", "Obesity", 
                              "Hypertensive diseases", "Thyroid gland disorders", 
                              "Other metabolic disorders", "Disorders of other endocrine glands") ~ "Cardiovascular, Vascular, and Metabolic Diseases",
    Full_Comorbidities %in% c("Metastatic cancer", "Other neoplasms") ~ "Metastatic Cancer and Neoplastic Diseases",
    TRUE ~ "Others" 
  ))

results_1 <- lapply(split(Data_1_indiv$Comorbidities_cato, f=Data_1_indiv$Baseline_medication), function(x){
  tmp <- table(x)
  tmp <- tmp[order(-tmp)]  # Sort by descending order of counts
  apply(rbind(tmp, round(tmp/sum(tmp)*100, digits=2)), 2, function(y){paste0(y[1], " (", y[2], "%)")})  # Format output
})


### Recoding 'Ethnicity' Variable
overall_clean_data$Ethnicity_recode <- recode(overall_clean_data$Ethnicity,
                                              "British" = "White",
                                              "Irish" = "White",
                                              "Any other white background" = "White",
                                              "White" = "White",
                                              .default = "Others")  # Recode 'Ethnicity' into broader categories

### Recoding 'Alcohol Intake Frequency'
overall_clean_data$Alcohol_intake_recode <- recode(overall_clean_data$Alcohol_intake_frequency,
                                                   "Never" = "Never",
                                                   "Special occasions only" = "Special occasions or 1-3 times a month",
                                                   "One to three times a month" = "Special occasions or 1-3 times a month",
                                                   "Once or twice a week" = "1-4 times a week",
                                                   "Three or four times a week" = "1-4 times a week",
                                                   "Daily or almost daily" = "Daily or almost daily")  

### Recoding and Merging IMD Values (Scotland, England, Wales)
Data_1_indiv <- Data_1_indiv %>%
  mutate(IMD_Scotland_Recode = case_when(
    IMD_Scotland < 18 ~ 1,
    IMD_Scotland >= 18 & IMD_Scotland <= 35.99 ~ 2,
    IMD_Scotland >= 36 & IMD_Scotland <= 53.99 ~ 3,
    IMD_Scotland >= 54 & IMD_Scotland <= 72 ~ 4,
    IMD_Scotland > 72 ~ 5
  ),
  IMD_England_Recode = case_when(
    IMD_England < 17 ~ 1,
    IMD_England >= 17 & IMD_England <= 33.99 ~ 2,
    IMD_England >= 34 & IMD_England <= 50.99 ~ 3,
    IMD_England >= 51 & IMD_England <= 67.99 ~ 4,
    IMD_England > 68 ~ 5
  ),
  IMD_Wales_Recode = case_when(
    IMD_Wales < 11 ~ 1,
    IMD_Wales >= 11 & IMD_Wales <= 21.99 ~ 2,
    IMD_Wales >= 22 & IMD_Wales <= 32.99 ~ 3,
    IMD_Wales >= 33 & IMD_Wales <= 43.99 ~ 4,
    IMD_Wales > 44 ~ 5
  ),
  Combined_IMD_Value = coalesce(IMD_England_Recode, IMD_Scotland_Recode, IMD_Wales_Recode))  # Merge IMD recoded values into a single column



### Filtering and Recoding 'Response_Type'
overall_clean_data <- overall_clean_data %>%
  filter(!is.na(Response_Type) & Response_Type != "Side effect" & Response_Type != "Unknown")  # Filter out NA, 'Side effect', and 'Unknown' values

# Recoding 'Response_Type' based on 'Baseline_medication' and 'End_of_study_medication'

overall_clean_data <- overall_clean_data %>%
  mutate(Response_Type = case_when(
    Baseline_medication == "Untreated" ~ "Untreated",
    Baseline_medication == "ACE inhibitors/ARBs" & End_of_study_medication == "ACE inhibitors/ARBs" ~ "ACE/ARBs_Responder",
    Baseline_medication == "ACE inhibitors/ARBs" & End_of_study_medication == "ACE inhibitors/ARBs + Calcium channel blockers" ~ "ACE/ARBs_Non_responder",
    Baseline_medication == "ACE inhibitors/ARBs" & End_of_study_medication == "ACE inhibitors/ARBs + Diuretics" ~ "ACE/ARBs_Non_responder",
    Baseline_medication == "ACE inhibitors/ARBs" & End_of_study_medication == "Untreated" ~ "Unknown",
    Baseline_medication == "ACE inhibitors/ARBs" & End_of_study_medication == "Calcium channel blockers" ~ "Side effect",
    Baseline_medication == "ACE inhibitors/ARBs" & is.na(End_of_study_medication) ~ NA_character_,
    Baseline_medication == "Calcium channel blockers" & End_of_study_medication == "Calcium channel blockers" ~ "CCB_Responder",
    Baseline_medication == "Calcium channel blockers" & End_of_study_medication == "Untreated" ~ "Unknown",
    Baseline_medication == "Calcium channel blockers" & is.na(End_of_study_medication) ~ NA_character_,
    Baseline_medication == "Calcium channel blockers" & End_of_study_medication == "ACE inhibitors/ARBs" ~ "Side effect",
    Baseline_medication == "Calcium channel blockers" & End_of_study_medication == "ACE inhibitors/ARBs + Calcium channel blockers" ~ "CCB_Non_responder",
    Baseline_medication == "Calcium channel blockers" & End_of_study_medication == "Calcium channel blockers + Diuretics" ~ "CCB_Non_responder",
    TRUE ~ "others"  # Optional: Captures any remaining cases not covered above
  ))  

# Convert 'Response_Type' back to a factor
overall_clean_data$Response_Type <- factor(overall_clean_data$Response_Type)  

# Count the number of missing values in 'Response_Type'
na_count <- sum(is.na(overall_clean_data$Response_Type))  
print(na_count)

### Saving Cleaned Data
write.csv(overall_clean_data, "/rds/general/project/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Full_clean_updated.csv", row.names = FALSE)  # Save the cleaned data
write.csv(overall_clean_updated, "/rds/general/project/hda_23-24/live/TDS/group3/TDS_group3/Dataset/Aim1_subset_data/overall_clean_data.csv", row.names = FALSE)  # Save the subset of cleaned data

### Adjusting and Cleaning 'Baseline_medication' Levels
Data_1_indiv <- subset(Data_1_indiv, Baseline_medication != "ACE inhibitors/ARBs + Calcium channel blockers + Diuretics")  # Remove specific levels
Data_1_indiv$Baseline_medication <- droplevels(Data_1_indiv$Baseline_medication)  

### Recoding 'Baseline_medication' for Simplification
Data_1_indiv$Baseline_medication_recode <- recode(Data_1_indiv$Baseline_medication,
                                                  "ACE inhibitors/ARBs" = "Step 1 (ACE/ARBs)",
                                                  "ACE inhibitors/ARBs + Calcium channel blockers" = "Step 2 (ACE/ARBs+CCB)",
                                                  "ACE inhibitors/ARBs + Diuretics" = "Step 2 (ACE/ARBs+Diuretics)",
                                                  "Calcium channel blockers" = "Step 1 (CCB)",
                                                  "Calcium channel blockers + Diuretics" = "Step 2 (CCB+Diuretics)",
                                                  "Untreated" = "Untreated")  # Simplify and recode 'Baseline_medication' levels

### Counting NA in 'End_of_study_medication'
na_count2 <- sum(is.na(Table_1B$End_of_study_medication))  
print(na_count2)
