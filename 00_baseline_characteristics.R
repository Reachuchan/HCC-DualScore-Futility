# ==============================================================================
# Project: HCC Futility Analysis - Surgical Resection vs. LRST Comparator
# Script Name: 01_table1_and_supp_s9_baseline_characteristics.R
# Author: Medical Bioinformatics & Clinical Research Group
# Target Output: 
#   - Table1_Baseline_Characteristics.xlsx
#   - Supplementary_Table_S9_Missing_Data_Summary.xlsx
# Description:
#   1. Raw data import, deduplication, and column name standardization.
#   2. Feature engineering (APRI, FIB-4, ALBI, DEI, OIDS, DeRitis Ratio, etc.).
#   3. Data cleaning and variable preparation using safe binary/factor mapping.
#   4. Construction and export of Table 1 via gtsummary.
#   5. Reverse Kaplan-Meier median follow-up estimation and manuscript draft text.
#   6. Full-dataset missing character inspection and Supplementary Table S9 export.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup & Library Loading
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(survival)
library(glmnet)
library(survminer)
library(mice)
library(Hmisc)
library(gtsummary)
library(writexl)


# ------------------------------------------------------------------------------
# 2. Resection Cohort: Data Import & Feature Engineering
# ------------------------------------------------------------------------------
df <- read_excel("HCC_futility _resection_analysis.xlsx")

# Deduplicate based on clinical columns excluding serial number ('sn')
df <- df %>% distinct(across(all_of(setdiff(names(df), "sn"))), .keep_all = TRUE)

# Standardize column names
safe_names <- colnames(df)
safe_names <- gsub(">=", "_gte_", safe_names)
safe_names <- gsub("\\+", "_plus_", safe_names)
safe_names <- gsub("\\.", "_", safe_names)
safe_names <- gsub(" ", "_", safe_names)
colnames(df) <- safe_names

# Clean pseudo-missing values and coerce continuous variables
df <- df %>%
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

# Define primary and secondary outcome endpoints
df <- df %>%
  mutate(
    O_time  = RFS_time_5y,
    O_event = RFS_status_5y,
    S_time  = ifelse(RFS_status_5y == 1, RFS_time_5y, OS_time_5y),
    S_event = ifelse(RFS_status_5y == 1, 0L, OS_status_5y)
  )

# ------------------------------------------------------------------------------
# 3. LRST Comparator Cohort: Data Import & Feature Engineering
# ------------------------------------------------------------------------------
df_LRST <- read_excel("HCC_futility_therapy_analysis.xlsx")

df_LRST <- df_LRST %>% distinct(across(all_of(setdiff(names(df_LRST), "sn"))), .keep_all = TRUE)

safe_names <- colnames(df_LRST)
safe_names <- gsub(">=", "_gte_", safe_names)
safe_names <- gsub("\\+", "_plus_", safe_names)
safe_names <- gsub("\\.", "_", safe_names)
safe_names <- gsub(" ", "_", safe_names)
colnames(df_LRST) <- safe_names

df_LRST <- df_LRST %>%
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

# ------------------------------------------------------------------------------
# 4. Data Cleaning & Variable Standardization Functions
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

df_clean      <- clean_and_coerce(df)
df_LRST_clean <- clean_and_coerce(df_LRST)

safe_binary <- function(x) {
  if (is.null(x)) return(rep(NA_real_, length(x)))
  x_str <- as.character(x)
  x_str <- gsub("\u00A0", "", x_str, fixed = TRUE)
  x_str <- trimws(x_str)
  case_when(
    is.na(x) ~ NA_real_,
    x_str %in% c(".", "..", '"."', "NA", "na", "NaN", "null", "NULL", "", " ", "unknown", "Unknown", "UNKNOWN") ~ NA_real_,
    x_str %in% c("1", "1.0", "TRUE", "true", "Y", "yes") ~ 1,
    x_str %in% c("0", "0.0", "FALSE", "false", "N", "no") ~ 0,
    TRUE ~ NA_real_
  )
}

prepare_vars <- function(data) {
  data %>%
    mutate(
      Age_num   = suppressWarnings(as.numeric(age)),
      AFP_num   = suppressWarnings(as.numeric(AFP)),
      PLT_num   = suppressWarnings(as.numeric(PLT)),
      AST_num   = suppressWarnings(as.numeric(AST)),
      ALT_num   = suppressWarnings(as.numeric(ALT)),
      Alb_num   = suppressWarnings(as.numeric(Alb)),
      Bili_num  = suppressWarnings(as.numeric(Bili)),
      ALBI_num  = suppressWarnings(as.numeric(ALBI_score)),
      Tumor_Dim = suppressWarnings(as.numeric(Diameter_of_largest_tumour_nodule)),
      
      Male             = safe_binary(gender),
      HBV              = safe_binary(hbv),
      HCV              = safe_binary(hcv),
      Cirrhosis        = safe_binary(cirr),
      Diabetes         = safe_binary(DM),
      Hypertension     = safe_binary(HTN),
      Coronary_disease = safe_binary(CAD),
      Chronic_kidney   = safe_binary(CKD),
      Ascites_var      = safe_binary(Ascites),
      
      AFP_gt_400       = ifelse(is.na(AFP_num), NA_real_, ifelse(AFP_num > 400, 1, 0)),
      Thrombocytopenia = ifelse(is.na(PLT_num), NA_real_, ifelse(PLT_num < 100, 1, 0)),
      Multinodular     = ifelse(is.na(No_of_Tumour_3), NA_real_, ifelse(No_of_Tumour_3 >= 2, 1, 0)),
      Rad_VI           = ifelse(is.na(Rad1VI), NA_real_, ifelse(Rad1VI == 1, 1, 0)),
      
      # Consolidate BCLC stage into a single factor variable
      BCLC_str   = toupper(trimws(as.character(BCLCStage))),
      BCLC_str   = gsub("\\.0$", "", BCLC_str),
      
      BCLC_Stage = case_when(
        BCLC_str %in% c("0", "STAGE 0", "VERY EARLY") ~ "Stage 0",
        BCLC_str %in% c("A", "1", "STAGE A")          ~ "Stage A",
        BCLC_str %in% c("B", "2", "B/C", "STAGE B")   ~ "Stage B",
        BCLC_str %in% c("C", "3", "STAGE C")          ~ "Stage C",
        BCLC_str %in% c("D", "4", "STAGE D")          ~ "Stage D",
        TRUE                                           ~ NA_character_
      ) %>% factor(levels = c("Stage 0", "Stage A", "Stage B", "Stage C", "Stage D")),
      
      Deaths   = ifelse(is.na(death_clean), NA_real_, ifelse(death_clean == 1, 1, 0))
    ) %>%
    select(-BCLC_str)
}

df_resection_prepared <- prepare_vars(df_clean) %>%
  mutate(Early_failure = ifelse(!is.na(O_event) & O_event == 1, 1, 0))

df_lrst_prepared <- prepare_vars(df_LRST_clean)

# ------------------------------------------------------------------------------
# 5. Build Table 1 (Baseline Characteristics) via gtsummary
# ------------------------------------------------------------------------------
tbl_resection <- df_resection_prepared %>%
  select(
    Age_num, Male, HBV, HCV, Cirrhosis, Diabetes, Ascites_var,
    Hypertension, Coronary_disease, Chronic_kidney,
    AFP_num, AFP_gt_400, PLT_num, Thrombocytopenia,
    AST_num, ALT_num, Alb_num, Bili_num, ALBI_num,
    Tumor_Dim, Multinodular, Rad_VI, BCLC_Stage, Deaths, Early_failure 
  ) %>%
  tbl_summary(
    missing = "no",
    percent = "column",
    type = list(
      c(Age_num, AFP_num, PLT_num, AST_num, ALT_num, Alb_num, Bili_num, ALBI_num, Tumor_Dim) ~ "continuous",
      c(Male, HBV, HCV, Cirrhosis, Diabetes, Hypertension, Coronary_disease, Chronic_kidney, Ascites_var, 
        AFP_gt_400, Thrombocytopenia, Multinodular, Rad_VI, Deaths, Early_failure) ~ "dichotomous",
      BCLC_Stage ~ "categorical"
    ),
    value = list(
      Male ~ 1, HBV ~ 1, HCV ~ 1, Cirrhosis ~ 1, Diabetes ~ 1, Ascites_var ~ 1,
      Hypertension ~ 1, Coronary_disease ~ 1, Chronic_kidney ~ 1,
      AFP_gt_400 ~ 1, Thrombocytopenia ~ 1, Multinodular ~ 1, Rad_VI ~ 1, 
      Deaths ~ 1, Early_failure ~ 1
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}–{p75})",
      all_dichotomous() ~ "{n} ({p}%)",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1, all_dichotomous() ~ c(0, 1), all_categorical() ~ c(0, 1)),
    label = list(
      Age_num ~ "Age, years (median, IQR)", Male ~ "Male sex, n (%)",
      HBV ~ "HBV, n (%)", HCV ~ "HCV, n (%)", Cirrhosis ~ "Cirrhosis, n (%)",
      Diabetes ~ "Diabetes, n (%)", Ascites_var ~ "Ascites, n (%)",
      Hypertension ~ "Hypertension, n (%)", Coronary_disease ~ "Coronary artery disease, n (%)",  
      Chronic_kidney ~ "Chronic kidney disease, n (%)", AFP_num ~ "AFP, ng/mL (median, IQR)",
      AFP_gt_400 ~ "AFP >400 ng/mL, n (%)", PLT_num ~ "Platelets ×10⁹/L (median, IQR)",
      Thrombocytopenia ~ "Thrombocytopenia (<100 ×10⁹/L), n (%)",
      AST_num ~ "AST, U/L (median, IQR)", ALT_num ~ "ALT, U/L (median, IQR)",
      Alb_num ~ "Albumin, g/L (median, IQR)", Bili_num ~ "Bilirubin, µmol/L (median, IQR)",
      ALBI_num ~ "ALBI score (median, IQR)", Tumor_Dim ~ "Largest tumour diameter, cm (median, IQR)",
      Multinodular ~ "Multinodular (≥2), n (%)", Rad_VI ~ "Radiologic vascular invasion, n (%)",
      BCLC_Stage ~ "BCLC stage, n (%)", Deaths ~ "Deaths during follow-up, n (%)", 
      Early_failure ~ "Early-failure events (resection only)*, n (%)"
    )
  )

tbl_lrst <- df_lrst_prepared %>%
  select(
    Age_num, Male, HBV, HCV, Cirrhosis, Diabetes, Ascites_var,
    Hypertension, Coronary_disease, Chronic_kidney,
    AFP_num, AFP_gt_400, PLT_num, Thrombocytopenia,
    AST_num, ALT_num, Alb_num, Bili_num, ALBI_num,
    Tumor_Dim, Multinodular, Rad_VI, BCLC_Stage, Deaths 
  ) %>%
  tbl_summary(
    missing = "no",
    percent = "column",
    type = list(
      c(Age_num, AFP_num, PLT_num, AST_num, ALT_num, Alb_num, Bili_num, ALBI_num, Tumor_Dim) ~ "continuous",
      c(Male, HBV, HCV, Cirrhosis, Diabetes, Hypertension, Coronary_disease, Chronic_kidney, Ascites_var, 
        AFP_gt_400, Thrombocytopenia, Multinodular, Rad_VI, Deaths) ~ "dichotomous",
      BCLC_Stage ~ "categorical"
    ),
    value = list(
      Male ~ 1, HBV ~ 1, HCV ~ 1, Cirrhosis ~ 1, Diabetes ~ 1, Ascites_var ~ 1,
      Hypertension ~ 1, Coronary_disease ~ 1, Chronic_kidney ~ 1,
      AFP_gt_400 ~ 1, Thrombocytopenia ~ 1, Multinodular ~ 1, Rad_VI ~ 1, 
      Deaths ~ 1
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}–{p75})",
      all_dichotomous() ~ "{n} ({p}%)",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1, all_dichotomous() ~ c(0, 1), all_categorical() ~ c(0, 1)),
    label = list(
      Age_num ~ "Age, years (median, IQR)", Male ~ "Male sex, n (%)",
      HBV ~ "HBV, n (%)", HCV ~ "HCV, n (%)", Cirrhosis ~ "Cirrhosis, n (%)",
      Diabetes ~ "Diabetes, n (%)", Ascites_var ~ "Ascites, n (%)",
      Hypertension ~ "Hypertension, n (%)", Coronary_disease ~ "Coronary artery disease, n (%)",  
      Chronic_kidney ~ "Chronic kidney disease, n (%)", AFP_num ~ "AFP, ng/mL (median, IQR)",
      AFP_gt_400 ~ "AFP >400 ng/mL, n (%)", PLT_num ~ "Platelets ×10⁹/L (median, IQR)",
      Thrombocytopenia ~ "Thrombocytopenia (<100 ×10⁹/L), n (%)",
      AST_num ~ "AST, U/L (median, IQR)", ALT_num ~ "ALT, U/L (median, IQR)",
      Alb_num ~ "Albumin, g/L (median, IQR)", Bili_num ~ "Bilirubin, µmol/L (median, IQR)",
      ALBI_num ~ "ALBI score (median, IQR)", Tumor_Dim ~ "Largest tumour diameter, cm (median, IQR)",
      Multinodular ~ "Multinodular (≥2), n (%)", Rad_VI ~ "Radiologic vascular invasion, n (%)",
      BCLC_Stage ~ "BCLC stage, n (%)", Deaths ~ "Deaths during follow-up, n (%)"
    )
  )

# Merge tables and format layout
table1_final <- tbl_merge(
  tbls = list(tbl_resection, tbl_lrst),
  tab_spanner = c("**Resection**", "**LRST comparator**")
) %>%
  bold_labels() %>%
  modify_footnote(
    everything() ~ "Continuous variables are presented as median (IQR) and categorical variables as n (%). Percentages were calculated based on available non-missing data (e.g., BCLC stage percentages use known stage cases as denominator). Complete distributions and counts of missing data for all baseline variables are provided in Supplementary Table S9."
  )

n_resection <- nrow(df_clean)
n_lrst      <- nrow(df_LRST_clean)

stat_cols <- grep("^stat_", colnames(table1_final$table_body), value = TRUE)

col_resection_title <- paste0("Resection (N=", n_resection, ")")
col_lrst_title      <- paste0("LRST comparator (N=", n_lrst, ")")

df_table1_output <- table1_final$table_body %>%
  select(label, all_of(stat_cols))

colnames(df_table1_output) <- c("Variable", col_resection_title, col_lrst_title)

df_table1_output <- df_table1_output %>%
  mutate(across(last_col(), ~ ifelse(is.na(.) | . == "", "—", .)))

n_row <- data.frame(
  Variable = "N",
  c2 = as.character(n_resection),
  c3 = as.character(n_lrst),
  stringsAsFactors = FALSE
)
colnames(n_row) <- colnames(df_table1_output)

df_table1_output <- bind_rows(n_row, df_table1_output)

# Export Table 1
write_xlsx(df_table1_output, "Table1_Baseline_Characteristics.xlsx")

# ------------------------------------------------------------------------------
# 6. Reverse Kaplan-Meier Follow-up Estimation & Results Draft
# ------------------------------------------------------------------------------
fit_fu_df <- survfit(
  Surv(as.numeric(OS_time_raw), as.numeric(death_clean) == 0) ~ 1,
  data = df_clean
)

fu_df <- summary(fit_fu_df)$table
cat(sprintf("Resection Median Follow-up: %.1f days (95%% CI: %.1f–%.1f days)\n", 
            fu_df["median"], fu_df["0.95LCL"], fu_df["0.95UCL"]))

fit_fu_lrst <- survfit(
  Surv(as.numeric(OS_time_raw), as.numeric(death_clean) == 0) ~ 1,
  data = df_LRST_clean
)

fu_lrst <- summary(fit_fu_lrst)$table
cat(sprintf("LRST Comparator Median Follow-up: %.1f days (95%% CI: %.1f–%.1f days)\n", 
            fu_lrst["median"], fu_lrst["0.95LCL"], fu_lrst["0.95UCL"]))

# Convert follow-up time to months and compute 5-year metrics
time_raw <- as.numeric(df_clean$OS_time_raw)
is_days  <- max(time_raw, na.rm = TRUE) > 100

if (is_days) {
  df_clean$OS_time_months <- time_raw / 30.4375
  t_5yr_val <- 60
} else {
  df_clean$OS_time_months <- time_raw
  t_5yr_val <- 5
}

df_clean$death_clean_num <- as.numeric(df_clean$death_clean)

fit_fu <- survfit(
  Surv(OS_time_months, death_clean_num == 0) ~ 1,
  data = df_clean
)
fu_sum <- summary(fit_fu)$table

fu_median <- round(fu_sum["median"], 1)
fu_lcl    <- round(fu_sum["0.95LCL"], 1)
fu_ucl    <- round(fu_sum["0.95UCL"], 1)

fit_os <- survfit(
  Surv(OS_time_months, death_clean_num) ~ 1,
  data = df_clean
)
os_5y_summary <- summary(fit_os, times = t_5yr_val)

n_at_risk_5y <- os_5y_summary$n.risk
total_n      <- nrow(df_clean)
pct_at_risk  <- round((n_at_risk_5y / total_n) * 100, 1)

total_deaths <- sum(df_clean$death_clean_num == 1, na.rm = TRUE)

manuscript_text <- sprintf(
  "The median follow-up duration was %.1f months (95%% CI, %.1f–%.1f months), and %s patients (%.1f%%) remained at risk at 5 years. A total of %s deaths from any cause were observed in the resection cohort.",
  fu_median, fu_lcl, fu_ucl,
  format(n_at_risk_5y, big.mark = ","),
  pct_at_risk,
  format(total_deaths, big.mark = ",")
)

cat("\n================================== Manuscript Results Draft ==================================\n\n")
cat(manuscript_text, "\n\n")
cat("===============================================================================================\n")

# ------------------------------------------------------------------------------
# 7. Supplementary Table S9 Construction & Export
# ------------------------------------------------------------------------------
get_rigorous_missing_summary <- function(data) {
  clean_data <- data %>%
    mutate(across(everything(), ~ {
      x_str <- as.character(.)
      x_str <- gsub("\u00A0", "", x_str, fixed = TRUE)
      x_str <- trimws(x_str)
      is_pseudo <- x_str %in% c(".", "..", '"."', "NA", "na", "NaN", "null", "NULL", "", " ", "unknown", "Unknown", "UNKNOWN")
      ifelse(is.na(.) | is_pseudo, NA_character_, x_str)
    }))
  
  missing_df <- data.frame(
    Variable           = colnames(clean_data),
    Missing_Count      = colSums(is.na(clean_data)),
    Missing_Percentage = round(colSums(is.na(clean_data)) / nrow(clean_data) * 100, 2),
    stringsAsFactors   = FALSE
  )
  missing_df[order(-missing_df$Missing_Count), ]
}

missing_summary_resection <- get_rigorous_missing_summary(df)
missing_summary_lrst      <- get_rigorous_missing_summary(df_LRST)

print(missing_summary_resection)
print(missing_summary_lrst)

clean_str <- function(x) {
  if (is.null(x)) return(rep(NA_character_, length(x)))
  x_str <- trimws(gsub("\u00A0", "", as.character(x), fixed = TRUE))
  is_pseudo <- x_str %in% c(".", "..", '"."', "NA", "na", "NaN", "null", "NULL", "", " ", "unknown", "Unknown", "UNKNOWN")
  ifelse(is.na(x) | is_pseudo, NA_character_, x_str)
}

get_missing <- function(data, col_name, is_albi = FALSE) {
  n_total <- nrow(data)
  if (is_albi) {
    vec_missing <- is.na(clean_str(data[["Alb"]])) | is.na(clean_str(data[["Bili"]]))
    m_n <- sum(vec_missing)
  } else {
    m_n <- sum(is.na(clean_str(data[[col_name]])))
  }
  sprintf("%s (%.1f%%)", format(m_n, big.mark = ","), round(m_n / n_total * 100, 1))
}

table1_vars <- list(
  list(label = "Age, years",                          col = "age"),
  list(label = "Male sex",                            col = "gender"),
  list(label = "Hepatitis B virus infection",         col = "hbv"),
  list(label = "Hepatitis C virus infection",         col = "hcv"),
  list(label = "Liver cirrhosis",                     col = "cirr"),
  list(label = "Diabetes mellitus",                   col = "DM"),
  list(label = "Hypertension",                        col = "HTN"),
  list(label = "Coronary artery disease",             col = "CAD"),
  list(label = "Chronic kidney disease",              col = "CKD"),
  list(label = "Ascites present",                     col = "Ascites"),
  list(label = "Serum AFP, ng/mL",                    col = "AFP"),
  list(label = "Platelet count, ×10⁹/L",              col = "PLT"),
  list(label = "Thrombocytopenia (<100 ×10⁹/L)",      col = "PLT"),
  list(label = "AST, U/L",                            col = "AST"),
  list(label = "ALT, U/L",                            col = "ALT"),
  list(label = "Albumin, g/L",                        col = "Alb"),
  list(label = "Total bilirubin, µmol/L",             col = "Bili"),
  list(label = "ALBI score",                          col = "Alb", is_albi = TRUE),
  list(label = "Largest tumour diameter, cm",         col_res = "Diameter of largest tumour nodule", col_thr = "Diameter_of_largest_tumour_nodule"),
  list(label = "Multinodular tumour (≥2 nodules)",    col_res = "No of Tumour 3",                    col_thr = "No_of_Tumour_3"),
  list(label = "Radiologic vascular invasion",        col = "Rad1VI"),
  list(label = "BCLC stage",                          col = "BCLCStage"),
  list(label = "All-cause mortality during follow-up",col = "death_clean")
)

supp9_df <- do.call(rbind, lapply(table1_vars, function(item) {
  c_res <- if (!is.null(item$col_res)) item$col_res else item$col
  c_thr <- if (!is.null(item$col_thr)) item$col_thr else item$col
  is_albi <- isTRUE(item$is_albi)
  
  data.frame(
    Variable = item$label,
    `Resection cohort (N=3,047) Missing, n (%)` = get_missing(df, c_res, is_albi),
    `LRST comparator cohort (N=711) Missing, n (%)` = get_missing(df_LRST, c_thr, is_albi),
    check.names = FALSE
  )
}))

# Export Supplementary Table S9
write_xlsx(supp9_df, "Supplementary_Table_S9_Missing_Data_Summary.xlsx")


