getwd()
wd = "/Users/sofi/Desktop/HDAICL/Term2/TDS3"
setwd(wd)
merge = load("merge_data.RData")


str(merge_data)
summary(merge_data)
table(merge_data$`HIV-1 seropositivity for Human Immunodeficiency Virus.0.0`)
summary(data$`Histology of cancer tumour.0.0`)

data = merge_data

colnames(data)

data$`month of birth.0.0`

names(data)[names(data) == "sex.0.0"] = "Sex"
names(data)[names(data) == "year of birth.0.0"] = "Year of birth"

#change column names
new_colnames = c(
  "EID",
  "Sex",
  "Year_of_birth",
  "Month_of_birth",
  "Date_of_attending_assessment_center",
  "Month_of_attending_assessment_centre",
  "Current_tobacco_smoking",
  "Past_tobacco_smoking",
  "Major_dietary_last_5_years",
  "Variation_in_diet",
  "Alcohol_intake_frequency",
  "Seen_doctor_for_nerves_anxiety_tension_or_depression",
  "Seen_psychiatrist_for_nerves_anxiety_tension_or_depression",
  "Had_major_operations",
  "Diabetes_diagnosed_by_doctor",
  "Cancer_diagnosed_by_doctor",
  "Number_of_depression_episodes",
  "Vascular_heart_problems_diagnosed_by_doctor",
  "Blood_clot_DVT_bronchitis_emphysema_asthma_rhinitis_eczema_allergy_diagnosed_by_doctor",
  "Smoking_status",
  "Ever_addicted_to_any_substance_or_behaviour",
  "Ever_addicted_to_alcohol",
  "Frequency_of_drinking_alcohol",
  "Ongoing_addiction_to_alcohol",
  "Ever_taken_cannabis",
  "Ever_addicted_to_illicit_or_recreational_drugs",
  "Ethnicity",
  "BMI",
  "Weight",
  "Age_when_attended_assessment_center",
  "Age_at_recruitment",
  "Doctor_diagnosed_COPD",
  "Age_COPD_diagnosed_by_doctor",
  "HIV-1_seropositivity_for_Human_Immunodeficiency_Virus",
  "Income_score_England",
  "Employment_score_England",
  "Health_score_England",
  "Education_score_England",
  "Housing_score_England",
  "Crime_score_England",
  "Living_environment_score_England",
  "Income_score_Wales",
  "Employment_score_Wales",
  "Health_score_Wales",
  "Education_score_Wales",
  "Access_to_services_score_Wales",
  "Housing_score_Wales",
  "Physical_environment_score_Wales",
  "Community_safety_score_Wales",
  "Index_of_Multiple_Deprivation_Wales",
  "Index_of_Multiple_Deprivation_Scotland",
  "Income_score_Scotland",
  "Health_score_Scotland",
  "Education_score_Scotland",
  "Housing_score_Scotland",
  "Access_to_services_score_Scotland",
  "Crime_score_Scotland",
  "Histology_of_cancer_tumour",
  "Type_of_cancer_ICD9",
  "Diagnoses_main_ICD10",
  "Diagnoses_main_ICD9",
  "Diagnoses_secondary_ICD10",
  "Diagnoses_secondary_ICD9",
  "Baseline_systolic_BP",
  "Baseline_diastolic_BP",
  "Baseline_medication",
  "Follow-up_BP_medication",
  "Follow-up_systolic_BP",
  "Follow-up_diastolic_BP",
  "End_of_study_medication",
  "Medication_status",
  "Hypertensive",
  "Controlled_baseline_BP"
)

colnames(data) = new_colnames

colnames(data)

#sum NAs per column
na_count = colSums(is.na(data))
na_count

table(data$Ever_addicted_to_any_substance_or_behaviour)
table(data$Ever_taken_cannabis)
table(data$Hypertensive)

#remove columns
data <- subset(data, select = -c(Month_of_attending_assessment_centre, Number_of_depression_episodes, 
                           Ever_addicted_to_alcohol, Frequency_of_drinking_alcohol, Ongoing_addiction_to_alcohol,
                           Age_when_attended_assessment_center))

table(data$Sex)
summary(data$Age_at_recruitment)
table(data$Type_of_cancer_ICD9)
table(data$Diagnoses_main_ICD10)

diag = as.data.frame(data$Diagnoses_main_ICD10)

getwd()
# Read the .tsv file
code10 = read.delim("diag10.tsv", header = TRUE, stringsAsFactors = FALSE)
code10 = subset(code10, select = -c(node_id, parent_id, selectable))

#add meanings to ICD codes--> diag10
diag10 = merge(data, code10, by.x = "Diagnoses_main_ICD10", by.y = "coding", all.x = TRUE)

table(diag10$ICD10_meaning)


#make EID first
diag10 = diag10[, c("EID", setdiff(names(diag10), "EID"))]
diag10$meaning

#change name for meaning to ICD10_meaning
colnames(diag10)[colnames(diag10) == "meaning"] <- "ICD10_meaning"
diag10$ICD10_meaning

#remove secondary diagnoses
diag10 = subset(diag10, select = -c(Diagnoses_secondary_ICD10, Diagnoses_secondary_ICD9))

#ICD9
code9 = read.delim("diag9.tsv", header = TRUE, stringsAsFactors = FALSE)
code9 = subset(code9, select = -c(node_id, parent_id, selectable))

#add meanings to ICD codes 9
diag109 = merge(diag10, code9, by.x = "Diagnoses_main_ICD9", by.y = "coding", all.x = TRUE)
colnames(diag109)[colnames(diag109) == "meaning"] <- "ICD9_meaning"

table(diag109$ICD9_meaning)

colnames(diag109)

df = subset(diag109, select = -c(Income_score_England, Employment_score_England, Health_score_England, Education_score_England, Housing_score_England, Crime_score_England, Living_environment_score_England, Income_score_Wales, Employment_score_Wales, Health_score_Wales, Education_score_Wales, Access_to_services_score_Wales, Housing_score_Wales, Physical_environment_score_Wales, Community_safety_score_Wales, Index_of_Multiple_Deprivation_Wales, Index_of_Multiple_Deprivation_Scotland, Income_score_Scotland, Health_score_Scotland, Education_score_Scotland, Housing_score_Scotland, Access_to_services_score_Scotland, Crime_score_Scotland))
colnames(df)


#add meanings to ICD codes--> diag10
diag10 = merge(data, code10, by.x = "Diagnoses_main_ICD10", by.y = "coding", all.x = TRUE)

table(diag10$ICD10_meaning)

#reorder columns
new_order = c("EID", setdiff(names(df), c("EID", "Diagnoses_main_ICD10", "Diagnoses_main_ICD9", "ICD10_meaning", "ICD9_meaning")), "Diagnoses_main_ICD10", "ICD10_meaning", "Diagnoses_main_ICD9", "ICD9_meaning")
final = df[new_order]

table(final$Sex)

write.csv(final, "final.csv", row.names = FALSE)
getwd()
str(final)

getwd()
wd = "/rds/general/project/hda_23-24/live/TDS/group3/TDS_group3/TDS_group3/Data"
setwd(wd)
getwd()
protein_merge <- readRDS("merge_data_with_protein.rds")


#IMD variable
getwd()
wd = "/rds/general/project/hda_23-24/live/TDS/group3/TDS_group3/extraction_and_recoding/outputs"
setwd(wd)
getwd()
my_data <- readRDS("ukb_extracted.rds")

my_data$`Index of Multiple Deprivation (Scotland).0.0`

new_data = my_data[, c("Index of Multiple Deprivation (Scotland).0.0", 
                        "Index of Multiple Deprivation (England).0.0", 
                        "Index of Multiple Deprivation (Wales).0.0")]



new_data$EID = rownames(new_data)

#merge index to final dataset
merge_final = merge(final, new_data, by = "EID", all = FALSE)

names(merge_final)[names(merge_final) == "Index of Multiple Deprivation (Scotland).0.0"] = "IMD_Scotland"
names(merge_final)[names(merge_final) == "Index of Multiple Deprivation (England).0.0"] = "IMD_England"
names(merge_final)[names(merge_final) == "Index of Multiple Deprivation (Wales).0.0"] = "IMD_Wales"


getwd()
wd = "/rds/general/project/hda_23-24/live/TDS/group3/TDS_group3/clean_dataset"
setwd(wd)
write.csv(merge_final, "merge_final.csv", row.names = FALSE)

summary(merge_final$IMD_Wales)
