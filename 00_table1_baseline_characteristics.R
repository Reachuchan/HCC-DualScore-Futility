# ==============================================================================
# Script Name: 00_table1_baseline_characteristics.R
# Purpose: Baseline characteristics table (Table 1) and median follow-up estimation 
#          for surgical resection and LRST comparator cohorts.
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
# 2. Resection Cohort: Data Import and Feature Engineering (Original Pipeline)
# ------------------------------------------------------------------------------
df <- read_excel("HCC_futility _resection_analysis.xlsx")
df <- df %>% distinct(across(all_of(setdiff(names(df), "sn"))), .keep_all = TRUE)

safe_names <- colnames(df)
safe_names <- gsub(">=", "_gte_", safe_names)
safe_names <- gsub("\\+", "_plus_", safe_names)
safe_names <- gsub("\\.", "_", safe_names)
safe_names <- gsub(" ", "_", safe_names)
colnames(df) <- safe_names

df_processed <- df %>%
  # Clean pseudo-NA characters
  mutate(across(everything(), ~ ifelse(. %in% c(".", "NaN", " ", "", "NA", "null", "NULL"), NA, .))) %>%
  
  # Coerce continuous endpoints and variables
  mutate(across(c(RFS_time_5y, RFS_status_5y, OS_time_5y, OS_status_5y,
                  age, BMI, PLT, INR, CR, Bili, Alb, AST, ALT, TBS, 
                  Diameter_of_largest_tumour_nodule, No_of_Tumour_Nodules, InitVascLymphInvasion), as.numeric)) %>%
  
  # Feature Engineering
  mutate(
    log_PLT        = log(PLT + 1),
    log_INR        = log(INR + 1),
    log_CR         = log(CR + 1),
    log_AST        = log(AST + 1),
    log_ALT        = log(ALT + 1),
    log_Bili       = log(Bili + 1),
    log_Alb        = log(Alb + 1),
    log_TBS        = log(TBS + 1),
    log_Diameter   = log(Diameter_of_largest_tumour_nodule + 1),
    log_Nodules    = log(No_of_Tumour_Nodules + 1),
    logAFP_plus_1     = as.numeric(`log(AFP_plus_1)`),
    No_of_Tumour_3    = as.numeric(as.character(No_of_Tumour_3)),
    FIB4_index        = (age * AST) / (PLT * sqrt(ALT)),
    log_FIB4          = log(FIB4_index + 1),
    log_CR            = log(CR + 1),
    log_INR           = log(INR + 1),
    ALBI_score        = (log10(pmax(Bili, 1e-6)) * 0.66) + (Alb * -0.085),
    Comorbidity_Count = (ifelse(CAD == 1, 1, 0) + ifelse(DM == 1, 1, 0) + ifelse(HTN == 1, 1, 0) + ifelse(CKD == 1, 1, 0)),
    Centre            = as.factor(centre_x)
  )

# Define Endpoints
df_resection <- df_processed %>%
  mutate(
    O_time  = RFS_time_5y,
    O_event = RFS_status_5y,
    S_time  = ifelse(RFS_status_5y == 1, RFS_time_5y, OS_time_5y),
    S_event = ifelse(RFS_status_5y == 1, 0L, OS_status_5y)
  )

# Prepare Resection Table 1 Display Variables
df_resection_prepared <- df_resection %>%
  mutate(
    Age_num          = age,
    AFP_num          = suppressWarnings(as.numeric(as.character(AFP))),
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
    
    Deaths           = ifelse(!is.na(death_clean) & death_clean == 1, 1, 0),
    Early_failure    = ifelse(!is.na(O_event) & O_event == 1, 1, 0)
  )

# ------------------------------------------------------------------------------
# 3. LRST Cohort: Data Import and Feature Engineering (Original Pipeline)
# ------------------------------------------------------------------------------
df_LRST <- read_excel("HCC_futility_therapy_analysis.xlsx")
df_LRST <- df_LRST %>% distinct(across(all_of(setdiff(names(df_LRST), "sn"))), .keep_all = TRUE)

safe_names_lrst <- colnames(df_LRST)
safe_names_lrst <- gsub(">=", "_gte_", safe_names_lrst)
safe_names_lrst <- gsub("\\+", "_plus_", safe_names_lrst)
safe_names_lrst <- gsub("\\.", "_", safe_names_lrst)
safe_names_lrst <- gsub(" ", "_", safe_names_lrst)
colnames(df_LRST) <- safe_names_lrst

df_lrst_processed <- df_LRST %>%
  mutate(across(everything(), ~ ifelse(. %in% c(".", "NaN", " ", "", "NA", "null", "NULL"), NA, .))) %>%
  mutate(across(c(age, PLT, gender, INR, CR, Bili, Alb, AST, ALT, 
                  Diameter_of_largest_tumour_nodule, OS_time_5y, OS_status_5y,
                  CAD, DM, HTN, CKD, Rad1VI, HBsAg, cirr, Ascites), ~ as.numeric(as.character(.)))) %>%
  mutate(
    log_Diameter      = log(Diameter_of_largest_tumour_nodule + 1),
    logAFP_plus_1     = as.numeric(`log(AFP_plus_1)`),
    No_of_Tumour_3    = as.numeric(as.character(No_of_Tumour_3)),
    FIB4_index        = (age * AST) / (PLT * sqrt(ALT)),
    log_FIB4          = log(FIB4_index + 1),
    log_CR            = log(CR + 1),
    log_INR           = log(INR + 1),
    ALBI_score        = (log10(pmax(Bili, 1e-6)) * 0.66) + (Alb * -0.085),
    Comorbidity_Count = (ifelse(CAD == 1, 1, 0) + ifelse(DM == 1, 1, 0) + ifelse(HTN == 1, 1, 0) + ifelse(CKD == 1, 1, 0)),
    Centre            = as.factor(centre_x)
  )

# Prepare LRST Table 1 Display Variables
df_lrst_prepared <- df_lrst_processed %>%
  mutate(
    Age_num          = age,
    AFP_num          = suppressWarnings(as.numeric(as.character(AFP))),
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

# ------------------------------------------------------------------------------
# 4. Generate and Export Table 1
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

# Export Excel
n_resection <- nrow(df_resection)
n_lrst      <- nrow(df_lrst_processed)

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
# 5. Reverse Kaplan-Meier Median Follow-up Estimation
# ------------------------------------------------------------------------------
if ("OS_time_raw" %in% colnames(df_resection) && "death_clean" %in% colnames(df_resection)) {
  fit_fu_df <- survfit(Surv(as.numeric(OS_time_raw), as.numeric(death_clean) == 0) ~ 1, data = df_resection)
  fu_df <- summary(fit_fu_df)$table
  cat(sprintf("Resection Median Follow-up: %.1f days (95%% CI: %.1f–%.1f days)\n", 
              fu_df["median"], fu_df["0.95LCL"], fu_df["0.95UCL"]))
}

if ("OS_time_raw" %in% colnames(df_lrst_processed) && "death_clean" %in% colnames(df_lrst_processed)) {
  fit_fu_lrst <- survfit(Surv(as.numeric(OS_time_raw), as.numeric(death_clean) == 0) ~ 1, data = df_lrst_processed)
  fu_lrst <- summary(fit_fu_lrst)$table
  cat(sprintf("LRST Comparator Median Follow-up: %.1f days (95%% CI: %.1f–%.1f days)\n", 
              fu_lrst["median"], fu_lrst["0.95LCL"], fu_lrst["0.95UCL"]))
}
