# ==============================================================================
# Script Name: 00_table1_baseline_characteristics.R
# Purpose: Clean and process baseline variables for the surgical resection cohort (n = 3,047)
#          and the LRST comparator cohort (n = 711). Generate Table 1 (gtsummary) and 
#          calculate reverse Kaplan-Meier median follow-up times for both cohorts.
# Outputs: Table1_Baseline_Characteristics.xlsx
# Language: R (>= 4.0.0)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(survival)
library(gtsummary)
library(writexl)

# ------------------------------------------------------------------------------
# 2. Helper Functions for Data Cleaning and Feature Coercion
# ------------------------------------------------------------------------------
clean_and_coerce <- function(data) {
  data %>%
    mutate(across(everything(), ~ ifelse(. %in% c(".", "NaN", " ", "", "NA", "null", "NULL"), NA_character_, as.character(.)))) %>%
    mutate(across(any_of(c(
      "age", "BMI", "PLT", "INR", "CR", "Bili", "Alb", "AST", "ALT", "AFP", "ALBI_score",
      "Diameter_of_largest_tumour_nodule", "No_of_Tumour_Nodules", "No_of_Tumour_3",
      "gender", "hbv", "hcv", "cirr", "DM", "HTN", "CAD", "CKD",
      "Ascites", "Rad1VI", "death_clean", "O_event"
    )), ~ suppressWarnings(as.numeric(.))))
}

prepare_vars <- function(data) {
  data %>%
    mutate(
      Age_num          = age,
      AFP_num          = AFP,
      PLT_num          = PLT,
      AST_num          = AST,
      ALT_num          = ALT,
      Alb_num          = Alb,
      Bili_num         = Bili,
      ALBI_num         = ALBI_score,
      Tumor_Dim        = Diameter_of_largest_tumour_nodule,
      
      Male             = ifelse(!is.na(gender) & gender == 1, 1, 0),
      HBV              = ifelse(!is.na(hbv) & hbv == 1, 1, 0),
      HCV              = ifelse(!is.na(hcv) & hcv == 1, 1, 0),
      Cirrhosis        = ifelse(!is.na(cirr) & cirr == 1, 1, 0),
      Diabetes         = ifelse(!is.na(DM) & DM == 1, 1, 0),
      Hypertension     = ifelse(!is.na(HTN) & HTN == 1, 1, 0),
      Coronary_disease = ifelse(!is.na(CAD) & CAD == 1, 1, 0),
      Chronic_kidney   = ifelse(!is.na(CKD) & CKD == 1, 1, 0),
      Ascites_var      = ifelse(!is.na(Ascites) & Ascites == 1, 1, 0),
      AFP_gt_400       = ifelse(!is.na(AFP_num) & AFP_num > 400, 1, 0),
      Thrombocytopenia = ifelse(!is.na(PLT_num) & PLT_num < 100, 1, 0),
      Multinodular     = ifelse(!is.na(No_of_Tumour_3) & No_of_Tumour_3 >= 2, 1, 0),
      Rad_VI           = ifelse(!is.na(Rad1VI) & Rad1VI == 1, 1, 0),
      BCLC_A           = ifelse(!is.na(BCLCStage) & toupper(as.character(BCLCStage)) %in% c("A", "1"), 1, 0),
      BCLC_B           = ifelse(!is.na(BCLCStage) & toupper(as.character(BCLCStage)) %in% c("B", "2"), 1, 0),
      BCLC_C           = ifelse(!is.na(BCLCStage) & toupper(as.character(BCLCStage)) %in% c("C", "3"), 1, 0),
      Deaths           = ifelse(!is.na(death_clean) & death_clean == 1, 1, 0)
    )
}

# ------------------------------------------------------------------------------
# 3. Load and Process Surgical Resection Cohort
# ------------------------------------------------------------------------------
df_resect_raw <- read_excel("HCC_futility_resection_analysis.xlsx")
df_resect_raw <- df_resect_raw %>% distinct(across(all_of(setdiff(names(df_resect_raw), "sn"))), .keep_all = TRUE)

safe_names <- colnames(df_resect_raw)
safe_names <- gsub(">=", "_gte_", safe_names); safe_names <- gsub("\\+", "_plus_", safe_names)
safe_names <- gsub("\\.", "_", safe_names); safe_names <- gsub(" ", "_", safe_names)
colnames(df_resect_raw) <- safe_names

df_resect_clean <- clean_and_coerce(df_resect_raw)
df_resection_prepared <- prepare_vars(df_resect_clean) %>%
  mutate(Early_failure = ifelse(!is.na(O_event) & O_event == 1, 1, 0))

# ------------------------------------------------------------------------------
# 4. Load and Process LRST Comparator Cohort
# ------------------------------------------------------------------------------
df_lrst_raw <- read_excel("HCC_futility_therapy_analysis.xlsx")
df_lrst_raw <- df_lrst_raw %>% distinct(across(all_of(setdiff(names(df_lrst_raw), "sn"))), .keep_all = TRUE)

safe_names_lrst <- colnames(df_lrst_raw)
safe_names_lrst <- gsub(">=", "_gte_", safe_names_lrst); safe_names_lrst <- gsub("\\+", "_plus_", safe_names_lrst)
safe_names_lrst <- gsub("\\.", "_", safe_names_lrst); safe_names_lrst <- gsub(" ", "_", safe_names_lrst)
colnames(df_lrst_raw) <- safe_names_lrst

df_lrst_clean <- clean_and_coerce(df_lrst_raw)
df_lrst_prepared <- prepare_vars(df_lrst_clean)

# ------------------------------------------------------------------------------
# 5. Generate Cohort-Specific gtsummary Tables
# ------------------------------------------------------------------------------
tbl_resection <- df_resection_prepared %>%
  select(
    Age_num, Male, HBV, HCV, Cirrhosis, Diabetes, Ascites_var,
    Hypertension, Coronary_disease, Chronic_kidney,
    AFP_num, AFP_gt_400, PLT_num, Thrombocytopenia,
    AST_num, ALT_num, Alb_num, Bili_num, ALBI_num,
    Tumor_Dim, Multinodular, Rad_VI, BCLC_A, BCLC_B, BCLC_C, Deaths, Early_failure 
  ) %>%
  tbl_summary(
    missing = "no",
    type = list(
      c(Age_num, AFP_num, PLT_num, AST_num, ALT_num, Alb_num, Bili_num, ALBI_num, Tumor_Dim) ~ "continuous",
      c(Male, HBV, HCV, Cirrhosis, Diabetes, Hypertension, Coronary_disease, Chronic_kidney, Ascites_var, AFP_gt_400, Thrombocytopenia, Multinodular, Rad_VI, BCLC_A, BCLC_B, BCLC_C, Deaths, Early_failure) ~ "dichotomous"
    ),
    value = list(
      Male ~ 1, HBV ~ 1, HCV ~ 1, Cirrhosis ~ 1, Diabetes ~ 1, Ascites_var ~ 1,
      Hypertension ~ 1, Coronary_disease ~ 1, Chronic_kidney ~ 1,
      AFP_gt_400 ~ 1, Thrombocytopenia ~ 1, Multinodular ~ 1, Rad_VI ~ 1,
      BCLC_A ~ 1, BCLC_B ~ 1, BCLC_C ~ 1, Deaths ~ 1, Early_failure ~ 1
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})", all_dichotomous() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_dichotomous() ~ c(0, 1)),
    label = list(
      Age_num ~ "Age, years (median, IQR)", Male ~ "Male sex, n (%)",
      HBV ~ "HBV, n (%)", HCV ~ "HCV, n (%)", Cirrhosis ~ "Cirrhosis, n (%)",
      Diabetes ~ "Diabetes, n (%)", Ascites_var ~ "Ascites, n (%)",
      Hypertension ~ "Hypertension, n (%)",                  
      Coronary_disease ~ "Coronary artery disease, n (%)",  
      Chronic_kidney ~ "Chronic kidney disease, n (%)",
      AFP_num ~ "AFP, ng/mL (median, IQR)", AFP_gt_400 ~ "AFP >400 ng/mL, n (%)",
      PLT_num ~ "Platelets ×10⁹/L (median, IQR)", Thrombocytopenia ~ "Thrombocytopenia (<100 ×10⁹/L), n (%)",
      AST_num ~ "AST, U/L (median, IQR)", ALT_num ~ "ALT, U/L (median, IQR)",
      Alb_num ~ "Albumin, g/L (median, IQR)", Bili_num ~ "Bilirubin, µmol/L (median, IQR)",
      ALBI_num ~ "ALBI score (median, IQR)", Tumor_Dim ~ "Largest tumour diameter, cm (median, IQR)",
      Multinodular ~ "Multinodular (≥2), n (%)", Rad_VI ~ "Radiologic vascular invasion, n (%)",
      BCLC_A ~ "BCLC stage A, n (%)", BCLC_B ~ "BCLC stage B, n (%)", BCLC_C ~ "BCLC stage C, n (%)",
      Deaths ~ "Deaths during follow-up, n (%)", Early_failure ~ "Early-failure events (resection only)*, n (%)"
    )
  )

tbl_lrst <- df_lrst_prepared %>%
  select(
    Age_num, Male, HBV, HCV, Cirrhosis, Diabetes, Ascites_var,
    Hypertension, Coronary_disease, Chronic_kidney,
    AFP_num, AFP_gt_400, PLT_num, Thrombocytopenia,
    AST_num, ALT_num, Alb_num, Bili_num, ALBI_num,
    Tumor_Dim, Multinodular, Rad_VI, BCLC_A, BCLC_B, BCLC_C, Deaths 
  ) %>%
  tbl_summary(
    missing = "no",
    type = list(
      c(Age_num, AFP_num, PLT_num, AST_num, ALT_num, Alb_num, Bili_num, ALBI_num, Tumor_Dim) ~ "continuous",
      c(Male, HBV, HCV, Cirrhosis, Diabetes, Hypertension, Coronary_disease, Chronic_kidney, Ascites_var, AFP_gt_400, Thrombocytopenia, Multinodular, Rad_VI, BCLC_A, BCLC_B, BCLC_C, Deaths) ~ "dichotomous"
    ),
    value = list(
      Male ~ 1, HBV ~ 1, HCV ~ 1, Cirrhosis ~ 1, Diabetes ~ 1, Ascites_var ~ 1,
      Hypertension ~ 1, Coronary_disease ~ 1, Chronic_kidney ~ 1,
      AFP_gt_400 ~ 1, Thrombocytopenia ~ 1, Multinodular ~ 1, Rad_VI ~ 1,
      BCLC_A ~ 1, BCLC_B ~ 1, BCLC_C ~ 1, Deaths ~ 1 
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})", all_dichotomous() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_dichotomous() ~ c(0, 1)),
    label = list(
      Age_num ~ "Age, years (median, IQR)", Male ~ "Male sex, n (%)",
      HBV ~ "HBV, n (%)", HCV ~ "HCV, n (%)", Cirrhosis ~ "Cirrhosis, n (%)",
      Diabetes ~ "Diabetes, n (%)", Ascites_var ~ "Ascites, n (%)",
      Hypertension ~ "Hypertension, n (%)",                  
      Coronary_disease ~ "Coronary artery disease, n (%)",  
      Chronic_kidney ~ "Chronic kidney disease, n (%)",
      AFP_num ~ "AFP, ng/mL (median, IQR)", AFP_gt_400 ~ "AFP >400 ng/mL, n (%)",
      PLT_num ~ "Platelets ×10⁹/L (median, IQR)", Thrombocytopenia ~ "Thrombocytopenia (<100 ×10⁹/L), n (%)",
      AST_num ~ "AST, U/L (median, IQR)", ALT_num ~ "ALT, U/L (median, IQR)",
      Alb_num ~ "Albumin, g/L (median, IQR)", Bili_num ~ "Bilirubin, µmol/L (median, IQR)",
      ALBI_num ~ "ALBI score (median, IQR)", Tumor_Dim ~ "Largest tumour diameter, cm (median, IQR)",
      Multinodular ~ "Multinodular (≥2), n (%)", Rad_VI ~ "Radiologic vascular invasion, n (%)",
      BCLC_A ~ "BCLC stage A, n (%)", BCLC_B ~ "BCLC stage B, n (%)", BCLC_C ~ "BCLC stage C, n (%)",
      Deaths ~ "Deaths during follow-up, n (%)"
    )
  )

# Merge Cohort Tables
table1_final <- tbl_merge(
  tbls = list(tbl_resection, tbl_lrst),
  tab_spanner = c("**Resection**", "**LRST comparator**")
) %>% bold_labels()

# Format and Export to Excel
n_resection <- nrow(df_resect_clean)
n_lrst      <- nrow(df_lrst_clean)

stat_cols <- grep("^stat_", colnames(table1_final$table_body), value = TRUE)
col_resection_title <- paste0("Resection (N=", n_resection, ")")
col_lrst_title      <- paste0("LRST comparator (N=", n_lrst, ")")

df_table1_output <- table1_final$table_body %>% select(label, all_of(stat_cols))
colnames(df_table1_output) <- c("Variable", col_resection_title, col_lrst_title)

df_table1_output <- df_table1_output %>% mutate(across(last_col(), ~ ifelse(is.na(.) | . == "", "—", .)))

n_row <- data.frame(
  Variable = "N",
  c2 = as.character(n_resection),
  c3 = as.character(n_lrst),
  stringsAsFactors = FALSE
)
colnames(n_row) <- colnames(df_table1_output)
df_table1_output <- bind_rows(n_row, df_table1_output)

write_xlsx(df_table1_output, "Table1_Baseline_Characteristics.xlsx")
cat(">> Table 1 successfully exported to Table1_Baseline_Characteristics.xlsx\n")

# ------------------------------------------------------------------------------
# 6. Reverse Kaplan-Meier Median Follow-up Estimation
# ------------------------------------------------------------------------------
fit_fu_df <- survfit(Surv(as.numeric(OS_time_raw), as.numeric(death_clean) == 0) ~ 1, data = df_resect_clean)
fu_df <- summary(fit_fu_df)$table
cat(sprintf("Resection Median Follow-up: %.1f days (95%% CI: %.1f–%.1f days)\n", 
            fu_df["median"], fu_df["0.95LCL"], fu_df["0.95UCL"]))

fit_fu_lrst <- survfit(Surv(as.numeric(OS_time_raw), as.numeric(death_clean) == 0) ~ 1, data = df_lrst_clean)
fu_lrst <- summary(fit_fu_lrst)$table
cat(sprintf("LRST Comparator Median Follow-up: %.1f days (95%% CI: %.1f–%.1f days)\n", 
            fu_lrst["median"], fu_lrst["0.95LCL"], fu_lrst["0.95UCL"]))
