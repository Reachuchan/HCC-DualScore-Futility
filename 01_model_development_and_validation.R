# ==============================================================================
# Script Name: 01_model_development_and_validation.R
# Purpose: Model training (Ridge Cox), bootstrap optimism correction, 
#          TRIPOD model cards, quadrant framework, LOCO cross-validation, 
#          and threshold grid sensitivity analyses.
# Outputs: Table 2, Table 4 (Panel A), Figure 1, Figure 3A, Figure 3D,
#          Supplementary Tables S1, S2, S7, S8, Supplementary Figures S5, S6.
# Language: R (>= 4.0.0)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(survival)
library(glmnet)
library(survminer)
library(mice)
library(Hmisc)
library(cmprsk)
library(patchwork)
library(grid)
library(ggExtra)
library(scales)
library(timeROC)
library(readr)

# ------------------------------------------------------------------------------
# 2. Data Import and Sanitation
# ------------------------------------------------------------------------------
df <- read_excel("HCC_futility _resection_analysis.xlsx")
df <- df %>% distinct(across(all_of(setdiff(names(df), "sn"))), .keep_all = TRUE)

safe_names <- colnames(df)
safe_names <- gsub(">=", "_gte_", safe_names)
safe_names <- gsub("\\+", "_plus_", safe_names)
safe_names <- gsub("\\.", "_", safe_names)
safe_names <- gsub(" ", "_", safe_names)
colnames(df) <- safe_names

# ------------------------------------------------------------------------------
# 3. Feature Engineering and Endpoint Definitions
# ------------------------------------------------------------------------------
df <- df %>%
  mutate(across(everything(), ~ ifelse(. == "." | . == "NaN" | . == " " | . == " ", NA, .))) %>%
  mutate(across(c(RFS_time_5y, RFS_status_5y, OS_time_5y, OS_status_5y,
                  age, BMI, PLT, INR, CR, Bili, Alb, AST, ALT, TBS, 
                  Diameter_of_largest_tumour_nodule, No_of_Tumour_Nodules, InitVascLymphInvasion), as.numeric)) %>%
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
    Comorbidity_Count = (ifelse(CAD == 1, 1, 0) + 
                         ifelse(DM == 1, 1, 0) + 
                         ifelse(HTN == 1, 1, 0) + 
                         ifelse(CKD == 1, 1, 0)),
    No_of_Tumour_3 = as.numeric(as.character(No_of_Tumour_3)),
    VascInv        = InitVascLymphInvasion,
    log_AST        = log(AST + 1),
    FIB4_index     = (age * AST) / (PLT * sqrt(ALT)),
    log_FIB4       = log(FIB4_index + 1),
    logAFP_plus_1  = `log(AFP_plus_1)`,
    ALBI_score     = (log10(pmax(Bili, 1e-6)) * 0.66) + (Alb * -0.085),
    Centre         = as.factor(centre_x)
  )

df <- df %>%
  mutate(
    O_time  = RFS_time_5y,
    O_event = RFS_status_5y,
    S_time  = ifelse(RFS_status_5y == 1, RFS_time_5y, OS_time_5y),
    S_event = ifelse(RFS_status_5y == 1, 0L, OS_status_5y)
  )

cat(sprintf(">> O-score Endpoint: %d events / %d total cases\n", sum(df$O_event, na.rm=TRUE), nrow(df)))
cat(sprintf(">> S-score Endpoint: %d events / %d total cases\n", sum(df$S_event, na.rm=TRUE), nrow(df)))

# ------------------------------------------------------------------------------
# 4. Predictor Definitions and MICE Imputation (m = 20)
# ------------------------------------------------------------------------------
O_VARS <- c("log_Diameter", "No_of_Tumour_3", "Rad1VI", "HBsAg","cirr","logAFP_plus_1")
S_VARS <- c("log_FIB4","log_CR", "log_INR","ALBI_score","Ascites","Comorbidity_Count")

cat(sprintf(">> Predictors defined: O-score (%d vars); S-score (%d vars)\n", length(O_VARS), length(S_VARS)))

all_vars <- unique(c(O_VARS, S_VARS,
                     "O_time","O_event","S_time","S_event",
                     "RFS_time_5y","RFS_status_5y",
                     "OS_time_5y","OS_status_5y",
                     "Centre"))

df_imp <- df %>% select(all_of(intersect(all_vars, colnames(df))))

cat("\n>> Missing Data Summary:\n")
miss_tbl <- data.frame(
  variable = names(df_imp),
  n_miss   = sapply(df_imp, function(x) sum(is.na(x))),
  pct_miss = round(sapply(df_imp, function(x) mean(is.na(x)))*100, 1)
) %>% filter(n_miss > 0) %>% arrange(desc(pct_miss))
print(miss_tbl, row.names=FALSE)

imp_method  <- make.method(df_imp)
outcome_cols <- c("O_time","O_event","S_time","S_event",
                  "RFS_time_5y","RFS_status_5y",
                  "OS_time_5y","OS_status_5y")
imp_method[intersect(outcome_cols, names(imp_method))] <- ""
imp_method["Centre"] <- ""

pred_matrix <- quickpred(df_imp, mincor=0.1, include=c("O_time","O_event","S_time","S_event"))
pred_matrix[intersect(c(outcome_cols,"Centre"), rownames(pred_matrix)), ] <- 0

cat(">> Running MICE (m=20, maxit=10)...\n")
imp_obj <- mice(df_imp, m=20, maxit=10, method=imp_method, predictorMatrix=pred_matrix, seed=2026, printFlag=FALSE)
cat(">> MICE Completed.\n")

dat_scored <- complete(imp_obj, 1)
dat_scored$O_time        <- df$O_time;        dat_scored$O_event  <- df$O_event
dat_scored$S_time        <- df$S_time;        dat_scored$S_event  <- df$S_event
dat_scored$OS_time_5y    <- df$OS_time_5y
dat_scored$OS_status_5y  <- df$OS_status_5y
dat_scored$RFS_time_5y   <- df$RFS_time_5y
dat_scored$RFS_status_5y <- df$RFS_status_5y
dat_scored$Centre        <- df$Centre

dat_scored <- dat_scored %>%
  filter(!is.na(O_time) & O_time > 0 & !is.na(O_event),
         !is.na(S_time) & S_time > 0 & !is.na(S_event))

cat(sprintf(">> Cleaned analytical cohort: %d cases\n", nrow(dat_scored)))

# ------------------------------------------------------------------------------
# 5. Core Fitting Helper Functions
# ------------------------------------------------------------------------------
get_lp_ridge <- function(train_dat, test_dat, x_vars, time_var, event_var, lambda_val = 0.5) {
  X_tr <- data.matrix(train_dat[, x_vars, drop = FALSE])
  X_te <- data.matrix(test_dat[, x_vars, drop = FALSE])
  y_tr <- Surv(train_dat[[time_var]], train_dat[[event_var]])
  
  for(j in seq_len(ncol(X_tr))) {
    m_val <- median(X_tr[, j], na.rm = TRUE)
    if(is.na(m_val)) m_val <- 0
    X_tr[is.na(X_tr[, j]), j] <- m_val
    X_te[is.na(X_te[, j]), j] <- m_val
  }
  
  ctr <- colMeans(X_tr)
  sds <- apply(X_tr, 2, sd)
  sds[sds == 0] <- 1
  
  X_tr_sc <- sweep(sweep(X_tr, 2, ctr, "-"), 2, sds, "/")
  X_te_sc <- sweep(sweep(X_te, 2, ctr, "-"), 2, sds, "/")
  
  fit <- glmnet(X_tr_sc, y_tr, family = "cox", alpha = 0, lambda = lambda_val, standardize = FALSE)
  raw_lp <- X_te_sc %*% coef(fit)
  lp_res <- unname(as.numeric(raw_lp))
  
  list(lp = lp_res, coefs = as.numeric(coef(fit)))
}

harrell_bootstrap <- function(data, x_vars, time_var, event_var, lambda_val = 0.5, B = 200, seed = 2026) {
  data <- as.data.frame(data)
  rownames(data) <- NULL
  n <- nrow(data)
  
  res_app <- get_lp_ridge(data, data, x_vars, time_var, event_var, lambda_val)
  lp_app  <- res_app$lp
  
  ok      <- !is.na(lp_app) & is.finite(lp_app) & data[[time_var]] > 0 & !is.na(data[[event_var]])
  time_v  <- unname(as.numeric(data[[time_var]][ok]))
  event_v <- unname(as.numeric(data[[event_var]][ok]))
  score_v <- unname(as.numeric(-lp_app[ok]))
  
  c_obj   <- concordance(Surv(time_v, event_v) ~ score_v)
  c_app   <- c_obj$concordance
  se_app  <- sqrt(c_obj$var)
  
  set.seed(seed)
  opt_vec <- corr_vec <- numeric(B)
  
  for (b in seq_len(B)) {
    idx <- sample(n, n, replace = TRUE)
    d_b <- data[idx, ]
    rownames(d_b) <- NULL
    
    lp_bb <- tryCatch(get_lp_ridge(d_b, d_b, x_vars, time_var, event_var, lambda_val)$lp, error = function(e) NULL)
    if (is.null(lp_bb)) { opt_vec[b] <- NA; corr_vec[b] <- NA; next }
    
    ok_bb    <- !is.na(lp_bb) & is.finite(lp_bb) & d_b[[time_var]] > 0 & !is.na(d_b[[event_var]])
    time_bb  <- unname(as.numeric(d_b[[time_var]][ok_bb]))
    event_bb <- unname(as.numeric(d_b[[event_var]][ok_bb]))
    score_bb <- unname(as.numeric(-lp_bb[ok_bb]))
    c_bb     <- concordance(Surv(time_bb, event_bb) ~ score_bb)$concordance
    
    lp_bo <- tryCatch(get_lp_ridge(d_b, data, x_vars, time_var, event_var, lambda_val)$lp, error = function(e) NULL)
    if (is.null(lp_bo)) { opt_vec[b] <- NA; corr_vec[b] <- NA; next }
    
    ok_bo    <- !is.na(lp_bo) & is.finite(lp_bo) & data[[time_var]] > 0 & !is.na(data[[event_var]])
    time_bo  <- unname(as.numeric(data[[time_var]][ok_bo]))
    event_bo <- unname(as.numeric(data[[event_var]][ok_bo]))
    score_bo <- unname(as.numeric(-lp_bo[ok_bo]))
    c_bo     <- concordance(Surv(time_bo, event_bo) ~ score_bo)$concordance
    
    opt_vec[b]  <- c_bb - c_bo
    corr_vec[b] <- c_app - (c_bb - c_bo)
  }
  
  opt_vec  <- opt_vec[!is.na(opt_vec)]
  corr_vec <- corr_vec[!is.na(corr_vec)]
  
  list(
    apparent  = c_app, app_lo = c_app - 1.96 * se_app, app_hi = c_app + 1.96 * se_app,
    optimism  = mean(opt_vec), corrected = c_app - mean(opt_vec),
    corr_lo   = quantile(corr_vec, 0.025), corr_hi = quantile(corr_vec, 0.975), sd = sd(corr_vec)
  )
}

cat("\n>> Running Bootstrap-200 Optimism Correction...\n")
res_o <- harrell_bootstrap(dat_scored, O_VARS, "O_time", "O_event")
res_s <- harrell_bootstrap(dat_scored, S_VARS, "S_time", "S_event")

metrics_summary <- data.frame(
  Model       = c(sprintf("O-score (%d-var)", length(O_VARS)), sprintf("S-score (%d-var)", length(S_VARS))),
  Apparent_C  = c(sprintf("%.3f (%.3f–%.3f)", res_o$apparent, res_o$app_lo, res_o$app_hi), sprintf("%.3f (%.3f–%.3f)", res_s$apparent, res_s$app_lo, res_s$app_hi)),
  Optimism    = c(sprintf("%.4f", res_o$optimism), sprintf("%.4f", res_s$optimism)),
  Corrected_C = c(sprintf("%.3f (%.3f–%.3f)", res_o$corrected, res_o$corr_lo, res_o$corr_hi), sprintf("%.3f (%.3f–%.3f)", res_s$corrected, res_s$corr_lo, res_s$corr_hi)),
  Boot_SD     = c(sprintf("%.4f", res_o$sd), sprintf("%.4f", res_s$sd))
)

cat("\n=== Metrics Summary Table ===\n")
print(metrics_summary, row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Table 2: Dual-Score Cox Model Coefficients
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Table 2]
var_labels <- c(
  "log_Diameter"      = "log(Tumour Diameter + 1)",
  "No_of_Tumour_3"    = "Tumour number (1, 2, ≥3)",
  "Rad1VI"            = "Radiologic vascular invasion",
  "HBsAg"             = "Chronic HBV (HBsAg positive)",
  "cirr"              = "Cirrhosis",
  "logAFP_plus_1"     = "log(AFP + 1)",
  "log_FIB4"          = "log(FIB-4 + 1)",
  "log_CR"            = "log(Serum creatinine + 1)",
  "log_INR"           = "log(INR + 1)",
  "ALBI_score"        = "ALBI score",
  "Ascites"           = "Ascites",
  "Comorbidity_Count" = "Cumulative comorbidity count (0–4)"
)

get_ridge_table_stats <- function(data, x_vars, time_var, event_var, lambda_val = 0.5, B = 200, seed = 2026) {
  X <- data.matrix(data[, x_vars, drop = FALSE])
  y <- Surv(data[[time_var]], data[[event_var]])
  
  ctr <- colMeans(X, na.rm = TRUE)
  sds <- apply(X, 2, sd, na.rm = TRUE); sds[sds == 0] <- 1
  X_sc <- sweep(sweep(X, 2, ctr, "-"), 2, sds, "/")
  
  fit_main <- glmnet(X_sc, y, family = "cox", alpha = 0, lambda = lambda_val, standardize = FALSE)
  beta_sc  <- as.numeric(coef(fit_main))
  beta_raw <- beta_sc / sds
  
  set.seed(seed)
  n <- nrow(data)
  boot_betas <- matrix(NA, nrow = B, ncol = length(x_vars))
  
  for (b in seq_len(B)) {
    idx <- sample(n, n, replace = TRUE)
    X_b <- X[idx, , drop = FALSE]
    y_b <- y[idx]
    
    ctr_b <- colMeans(X_b, na.rm = TRUE)
    sds_b <- apply(X_b, 2, sd, na.rm = TRUE); sds_b[sds_b == 0] <- 1
    X_b_sc <- sweep(sweep(X_b, 2, ctr_b, "-"), 2, sds_b, "/")
    
    fit_b <- tryCatch(glmnet(X_b_sc, y_b, family = "cox", alpha = 0, lambda = lambda_val, standardize = FALSE), error = function(e) NULL)
    if (!is.null(fit_b)) {
      boot_betas[b, ] <- as.numeric(coef(fit_b)) / sds_b
    }
  }
  
  se_vec <- apply(boot_betas, 2, sd, na.rm = TRUE)
  se_vec[se_vec == 0 | is.na(se_vec)] <- 0.001
  
  z_scores <- beta_raw / se_vec
  p_vals   <- 2 * (1 - pnorm(abs(z_scores)))
  hr_val   <- exp(beta_raw)
  hr_lo    <- exp(beta_raw - 1.96 * se_vec)
  hr_hi    <- exp(beta_raw + 1.96 * se_vec)
  
  res_df <- data.frame(
    Variable = x_vars,
    beta     = sprintf("%.3f", beta_raw),
    SE       = sprintf("%.3f", se_vec),
    HR       = sprintf("%.2f", hr_val),
    CI_95    = sprintf("%.2f–%.2f", hr_lo, hr_hi),
    p_val    = ifelse(p_vals < 0.001, "<0.001", sprintf("%.3f", p_vals)),
    stringsAsFactors = FALSE
  )
  res_df$Variable <- ifelse(res_df$Variable %in% names(var_labels), var_labels[res_df$Variable], res_df$Variable)
  n_events <- sum(data[[event_var]] == 1, na.rm = TRUE)
  
  list(df = res_df, n_events = n_events)
}

res_o_ridge <- get_ridge_table_stats(dat_scored, O_VARS, "O_time", "O_event", lambda_val = 0.5)
res_s_ridge <- get_ridge_table_stats(dat_scored, S_VARS, "S_time", "S_event", lambda_val = 0.5)

header_o <- data.frame(Variable = "O-score (recurrence)", beta = "", SE = "", HR = "", CI_95 = "", p_val = "", stringsAsFactors = FALSE)
header_s <- data.frame(Variable = "S-score (non-recurrent mortality)", beta = "", SE = "", HR = "", CI_95 = "", p_val = "", stringsAsFactors = FALSE)

Table2_Ridge <- bind_rows(header_o, res_o_ridge$df, header_s, res_s_ridge$df)
colnames(Table2_Ridge) <- c("Variable", "β", "SE", "HR", "95% CI", "p-value")

cat("\n=== Table 2. Dual-score Cox model coefficients ===\n")
print(Table2_Ridge, row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Baseline Hazard and Risk Score Calculations
# ------------------------------------------------------------------------------
Xo <- data.matrix(dat_scored[, O_VARS])
for(j in seq_len(ncol(Xo))) Xo[is.na(Xo[,j]),j] <- median(Xo[,j], na.rm=TRUE)
ctr_o <- colMeans(Xo); sd_o <- apply(Xo,2,sd); sd_o[sd_o==0] <- 1
Xo_sc <- sweep(sweep(Xo,2,ctr_o,"-"),2,sd_o,"/")
fit_o <- glmnet(Xo_sc, Surv(dat_scored$O_time, dat_scored$O_event), family="cox", alpha=0, lambda=0.5, standardize=FALSE)
dat_scored$O_lp <- as.numeric(Xo_sc %*% coef(fit_o))

Xs <- data.matrix(dat_scored[, S_VARS])
for(j in seq_len(ncol(Xs))) Xs[is.na(Xs[,j]),j] <- median(Xs[,j], na.rm=TRUE)
ctr_s <- colMeans(Xs); sd_s <- apply(Xs,2,sd); sd_s[sd_s==0] <- 1
Xs_sc <- sweep(sweep(Xs,2,ctr_s,"-"),2,sd_s,"/")
fit_s <- glmnet(Xs_sc, Surv(dat_scored$S_time, dat_scored$S_event), family="cox", alpha=0, lambda=0.5, standardize=FALSE)
dat_scored$S_lp <- as.numeric(Xs_sc %*% coef(fit_s))

t5yr <- 1826.25
bh_o <- basehaz(coxph(Surv(O_time,O_event)~offset(O_lp), data=dat_scored[dat_scored$O_time>0,], ties="efron"), centered=FALSE)
bh_s <- basehaz(coxph(Surv(S_time,S_event)~offset(S_lp), data=dat_scored[dat_scored$S_time>0,], ties="efron"), centered=FALSE)
h0_o <- max(bh_o$hazard[bh_o$time <= t5yr], na.rm=TRUE)
h0_s <- max(bh_s$hazard[bh_s$time <= t5yr], na.rm=TRUE)

dat_scored <- dat_scored %>%
  mutate(pred_O_5yr = exp(-h0_o * exp(O_lp)),
         pred_S_5yr = exp(-h0_s * exp(S_lp)))

coef_tbl_o <- data.frame(Variable = O_VARS, Coefficient = as.numeric(coef(fit_o)))
coef_tbl_s <- data.frame(Variable = S_VARS, Coefficient = as.numeric(coef(fit_s)))

# ------------------------------------------------------------------------------
# 8. Supplementary Tables S1 & S2: TRIPOD-Compliant Model Cards
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Tables S1 and S2]
coding_dict <- c(
  "log_Diameter"      = "Tumour diameter, log(cm + 1)",
  "No_of_Tumour_3"    = "Ordinal numeric coding: 1 = single tumour, 2 = two tumours, 3 = ≥3 tumours",
  "Rad1VI"            = "Binary indicator: radiologic vascular invasion (0 = No, 1 = Yes)",
  "HBsAg"             = "Binary indicator: chronic HBV status (0 = Negative, 1 = Positive)",
  "cirr"              = "Binary indicator: liver cirrhosis (0 = No, 1 = Yes)",
  "logAFP_plus_1"     = "Serum AFP level, log(ng/mL + 1)",
  "log_FIB4"          = "FIB-4 index, log(value + 1)",
  "log_CR"            = "Serum creatinine, log(μmol/L + 1)",
  "log_INR"           = "Prothrombin time INR, log(value + 1)",
  "ALBI_score"        = "ALBI score, continuous variable",
  "Ascites"           = "Binary indicator: clinical ascites (0 = None, 1 = Present)",
  "Comorbidity_Count" = "Count variable: cumulative comorbidities (0, 1, 2, ≥3)"
)

generate_tripod_model_card <- function(data, x_vars, time_var, event_var, fit_glmnet, lambda_val = 0.5, B = 200, seed = 2026) {
  X_raw <- data.matrix(data[, x_vars, drop = FALSE])
  y     <- Surv(data[[time_var]], data[[event_var]])
  n     <- nrow(data)
  
  sds <- apply(X_raw, 2, sd, na.rm = TRUE); sds[sds == 0] <- 1
  beta_sc  <- as.numeric(coef(fit_glmnet))
  beta_raw <- beta_sc / sds
  hr_point <- exp(beta_raw)
  
  set.seed(seed)
  boot_betas <- matrix(NA, nrow = B, ncol = length(x_vars))
  idx_events    <- which(data[[event_var]] == 1)
  idx_nonevents <- which(data[[event_var]] == 0)
  n_events      <- length(idx_events)
  n_nonevents   <- length(idx_nonevents)
  
  for (b in seq_len(B)) {
    boot_events    <- sample(idx_events, size = n_events, replace = TRUE)
    boot_nonevents <- sample(idx_nonevents, size = n_nonevents, replace = TRUE)
    idx            <- c(boot_events, boot_nonevents)
    
    X_b <- X_raw[idx, , drop = FALSE]; y_b <- y[idx]
    ctr_b <- colMeans(X_b, na.rm = TRUE)
    sds_b <- apply(X_b, 2, sd, na.rm = TRUE); sds_b[sds_b == 0] <- 1
    X_b_sc <- sweep(sweep(X_b, 2, ctr_b, "-"), 2, sds_b, "/")
    
    fit_b <- tryCatch(glmnet(X_b_sc, y_b, family = "cox", alpha = 0, lambda = lambda_val, standardize = FALSE), error = function(e) NULL)
    if (!is.null(fit_b)) boot_betas[b, ] <- as.numeric(coef(fit_b)) / sds_b
  }
  
  ci_lo_beta <- apply(boot_betas, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_hi_beta <- apply(boot_betas, 2, quantile, probs = 0.975, na.rm = TRUE)
  ci_lo_hr   <- exp(ci_lo_beta)
  ci_hi_hr   <- exp(ci_hi_beta)
  
  var_clean_names <- c(
    "log_Diameter" = "Tumour diameter", "No_of_Tumour_3" = "Tumour number", "Rad1VI" = "Radiologic vascular invasion",
    "HBsAg" = "Chronic HBV", "cirr" = "Cirrhosis", "logAFP_plus_1" = "Serum AFP", "log_FIB4" = "FIB-4 index",
    "log_CR" = "Serum creatinine", "log_INR" = "INR", "ALBI_score" = "ALBI score", "Ascites" = "Ascites",
    "Comorbidity_Count" = "Comorbidity count"
  )
  
  top_table <- data.frame(
    `Predictor`        = ifelse(x_vars %in% names(var_clean_names), var_clean_names[x_vars], x_vars),
    `Coding`           = ifelse(x_vars %in% names(coding_dict), coding_dict[x_vars], x_vars),
    `β`                = sprintf("%.4f", beta_raw),
    `HR`               = sprintf("%.2f", hr_point),
    `Bootstrap 95% CI` = sprintf("%.2f–%.2f", ci_lo_hr, ci_hi_hr),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  
  lp_raw <- as.numeric(X_raw %*% beta_raw)
  sub_data <- data; sub_data$lp_raw_calc <- lp_raw
  sub_data <- sub_data %>% filter(!is.na(.data[[time_var]]) & .data[[time_var]] > 0 & !is.na(.data[[event_var]]))
  
  fit_cox_offset <- coxph(Surv(sub_data[[time_var]], sub_data[[event_var]]) ~ offset(lp_raw_calc), data = sub_data, ties = "efron")
  bh <- basehaz(fit_cox_offset, centered = FALSE)
  
  max_t <- max(sub_data[[time_var]], na.rm = TRUE)
  target_days <- if (max_t > 100) c(730.5, 1095.75, 1826.25) else c(2.0, 3.0, 5.0)
  horizons    <- c("2y", "3y", "5y")
  h0_vals <- numeric(3); s0_vals <- numeric(3)
  
  for (i in seq_along(target_days)) {
    t_target <- target_days[i]
    valid_h  <- bh$hazard[bh$time <= t_target]
    h0_val   <- if (length(valid_h) > 0) max(valid_h, na.rm = TRUE) else 0
    h0_vals[i] <- h0_val
    s0_vals[i] <- exp(-h0_val)
  }
  
  bottom_table <- data.frame(
    `Horizon` = horizons, `H0(t)` = sprintf("%.6f", h0_vals), `S0(t)` = sprintf("%.6f", s0_vals),
    check.names = FALSE
  )
  list(top = top_table, bottom = bottom_table, n_total = nrow(sub_data), n_events = sum(sub_data[[event_var]] == 1, na.rm = TRUE))
}

card_o_final <- generate_tripod_model_card(dat_scored, O_VARS, "O_time", "O_event", fit_o, lambda_val = 0.5)
card_s_final <- generate_tripod_model_card(dat_scored, S_VARS, "S_time", "S_event", fit_s, lambda_val = 0.5)

cat("\n=== Supplementary Table S1: O-score Model Card ===\n")
print(card_o_final$top, row.names = FALSE)
print(card_o_final$bottom, row.names = FALSE)

cat("\n=== Supplementary Table S2: S-score Model Card ===\n")
print(card_s_final$top, row.names = FALSE)
print(card_s_final$bottom, row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. Quadrant Classification and Summary Analysis
# ------------------------------------------------------------------------------
O_pct_cutoff <- 0.30
S_pct_cutoff <- 0.15

O_cutoff_pct <- quantile(dat_scored$pred_O_5yr, O_pct_cutoff, na.rm=TRUE)
S_cutoff_pct <- quantile(dat_scored$pred_S_5yr, S_pct_cutoff, na.rm=TRUE)
O_cutoff_abs <- O_cutoff_pct
S_cutoff_abs <- S_cutoff_pct

dat_scored <- dat_scored %>%
  mutate(
    high_O = as.integer(pred_O_5yr < O_cutoff_abs),
    high_S = as.integer(pred_S_5yr < S_cutoff_abs),
    quadrant = case_when(
      high_O==0 & high_S==0 ~ "Low-O / Low-S",
      high_O==1 & high_S==0 ~ "High-O / Low-S",
      high_O==0 & high_S==1 ~ "Low-O / High-S",
      high_O==1 & high_S==1 ~ "High-O / High-S",
      TRUE ~ NA_character_
    ),
    quadrant = factor(quadrant, levels=c("Low-O / Low-S","High-O / Low-S","Low-O / High-S","High-O / High-S"))
  )

target_cause_O <- 1
target_cause_S <- if (2 %in% unique(dat_scored$S_event)) 2 else 1

get_cif_at_5yr <- function(time_vec, event_vec, target_cause, t_target) {
  if (!any(event_vec == target_cause, na.rm = TRUE)) return(0)
  ci <- tryCatch(cuminc(ftime = time_vec, fstatus = event_vec, cencode = 0), error = function(e) NULL)
  if (is.null(ci)) return(0)
  curve_keys <- setdiff(names(ci), "Tests")
  if (length(curve_keys) == 0) return(0)
  
  matching_key <- NULL
  for (k in curve_keys) {
    tokens <- unlist(strsplit(k, " "))
    if (tokens[length(tokens)] == as.character(target_cause)) { matching_key <- k; break }
  }
  if (is.null(matching_key) && length(curve_keys) == 1) matching_key <- curve_keys[1]
  if (is.null(matching_key)) return(0)
  
  t_vals <- ci[[matching_key]]$time
  est_vals <- ci[[matching_key]]$est
  idx <- max(which(t_vals <= t_target))
  if (is.finite(idx) && idx > 0) return(est_vals[idx] * 100) else return(0)
}

quad_summary_final <- data.frame()
t_5yr <- ifelse(max(c(dat_scored$O_time, dat_scored$S_time, dat_scored$OS_time_5y), na.rm = TRUE) > 100, 1826.25, 5)

for (q in levels(dat_scored$quadrant)) {
  d_q <- dat_scored %>% filter(quadrant == q)
  n_q <- nrow(d_q)
  pct_q <- n_q / nrow(dat_scored) * 100
  
  fit_os <- survfit(Surv(OS_time_5y, OS_status_5y) ~ 1, data = d_q)
  sum_os <- summary(fit_os, times = t_5yr)
  os_5yr <- if (length(sum_os$surv) > 0) sum_os$surv * 100 else NA
  os_lower <- if (length(sum_os$lower) > 0) sum_os$lower * 100 else NA
  os_upper <- if (length(sum_os$upper) > 0) sum_os$upper * 100 else NA
  
  cif_rec <- get_cif_at_5yr(d_q$O_time, d_q$O_event, target_cause = target_cause_O, t_target = t_5yr)
  cif_nonrec <- get_cif_at_5yr(d_q$S_time, d_q$S_event, target_cause = target_cause_S, t_target = t_5yr)
  
  quad_summary_final <- rbind(quad_summary_final, data.frame(
    `Quadrant Name` = q, `Count_N` = n_q, `Pct` = sprintf("%.1f%%", pct_q),
    `5-yr OS (95% CI)` = sprintf("%.1f%% (%.1f%%-%.1f%%)", os_5yr, os_lower, os_upper),
    `5-yr Recurrence CIF` = sprintf("%.1f%%", cif_rec),
    `5-yr Non-recurrent Death CIF` = sprintf("%.1f%%", cif_nonrec),
    check.names = FALSE
  ))
}

# ------------------------------------------------------------------------------
# 10. Figure 1: Dual-Score Quadrant Framework (2x2 Display)
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Figure 1]
get_plot_data <- function(df, q_name) {
  row <- df[df$`Quadrant Name` == q_name, ]
  n_str <- sprintf("n = %s (%s)", format(row$Count_N, big.mark = ","), row$Pct)
  raw_os_ci <- row$`5-yr OS (95% CI)`
  parts <- strsplit(raw_os_ci, " ")[[1]]
  os_main <- parts[1]
  os_ci_str <- paste0("(95% CI: ", gsub("^\\(|\\)$", "", parts[2]), ")")
  list(n = n_str, os = os_main, ci = os_ci_str, rec = row$`5-yr Recurrence CIF`, nonrec = row$`5-yr Non-recurrent Death CIF`)
}

ll <- get_plot_data(quad_summary_final, "Low-O / Low-S")
hl <- get_plot_data(quad_summary_final, "High-O / Low-S")
lh <- get_plot_data(quad_summary_final, "Low-O / High-S")
hh <- get_plot_data(quad_summary_final, "High-O / High-S")

rect_df <- data.frame(
  xmin = c(0.08, 1.02, 0.08, 1.02), xmax = c(0.98, 1.92, 0.98, 1.92),
  ymin = c(0.06, 0.06, 1.02, 1.02), ymax = c(0.98, 0.98, 1.94, 1.94),
  fill_color = c("#A6D9AC", "#A6D9AC", "#FCD281", "#F7A5A5")
)

x_ctr_col1 <- (0.08 + 0.98) / 2; x_left_col1 <- 0.11
x_ctr_col2 <- (1.02 + 1.92) / 2; x_left_col2 <- 1.05

p_fig1 <- ggplot() +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color), color = "black", linewidth = 1.2) +
  scale_fill_identity() +
  
  annotate("text", x = x_ctr_col1, y = 0.91, label = "Low O / Low S", fontface = "bold", size = 5.0, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 0.85, label = "(standard resection)", fontface = "italic", size = 3.7, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 0.73, label = ll$n, size = 3.9, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 0.61, label = "5-yr OS:", size = 3.6, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 0.49, label = ll$os, fontface = "bold", size = 6, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 0.38, label = ll$ci, fontface = "italic", size = 3.2, hjust = 0.5) +
  annotate("text", x = x_left_col1, y = 0.27, label = paste0("5-yr CIF recurrence: ", ll$rec), size = 3.3, hjust = 0) +
  annotate("text", x = x_left_col1, y = 0.20, label = paste0("5-yr CIF non-rec death: ", ll$nonrec), size = 3.3, hjust = 0) +
  annotate("text", x = x_left_col1, y = 0.13, label = "ATT 5-yr RD vs LRST: +20.2 (+7.3, +34.3) pp", fontface = "italic", size = 3.1, hjust = 0) +
  
  annotate("text", x = x_ctr_col2, y = 0.91, label = "High O / Low S", fontface = "bold", size = 5.0, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 0.85, label = "(oncological futility)", fontface = "italic", size = 3.7, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 0.73, label = hl$n, size = 3.9, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 0.61, label = "5-yr OS:", size = 3.6, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 0.49, label = hl$os, fontface = "bold", size = 6, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 0.38, label = hl$ci, fontface = "italic", size = 3.2, hjust = 0.5) +
  annotate("text", x = x_left_col2, y = 0.27, label = paste0("5-yr CIF recurrence: ", hl$rec), size = 3.3, hjust = 0) +
  annotate("text", x = x_left_col2, y = 0.20, label = paste0("5-yr CIF non-rec death: ", hl$nonrec), size = 3.3, hjust = 0) +
  annotate("text", x = x_left_col2, y = 0.13, label = "ATT 5-yr RD vs LRST: +6.4 (-5.4, +20.7) pp", fontface = "italic", size = 3.1, hjust = 0) +
  
  annotate("text", x = x_ctr_col1, y = 1.87, label = "Low O / High S", fontface = "bold", size = 5.0, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 1.81, label = "(surgical futility)", fontface = "italic", size = 3.7, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 1.69, label = lh$n, size = 3.9, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 1.57, label = "5-yr OS:", size = 3.6, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 1.45, label = lh$os, fontface = "bold", size = 6, hjust = 0.5) +
  annotate("text", x = x_ctr_col1, y = 1.34, label = lh$ci, fontface = "italic", size = 3.2, hjust = 0.5) +
  annotate("text", x = x_left_col1, y = 1.23, label = paste0("5-yr CIF recurrence: ", lh$rec), size = 3.3, hjust = 0) +
  annotate("text", x = x_left_col1, y = 1.16, label = paste0("5-yr CIF non-rec death: ", lh$nonrec), size = 3.3, hjust = 0) +
  annotate("text", x = x_left_col1, y = 1.09, label = "ATT 5-yr RD vs LRST: +18.8 (+3.7, +34.2) pp", fontface = "italic", size = 3.1, hjust = 0) +
  
  annotate("text", x = x_ctr_col2, y = 1.87, label = "High O / High S", fontface = "bold", size = 5.0, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 1.81, label = "(doomed)", fontface = "italic", size = 3.7, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 1.69, label = hh$n, size = 3.9, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 1.57, label = "5-yr OS:", size = 3.6, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 1.45, label = hh$os, fontface = "bold", size = 6, hjust = 0.5) +
  annotate("text", x = x_ctr_col2, y = 1.34, label = hh$ci, fontface = "italic", size = 3.2, hjust = 0.5) +
  annotate("text", x = x_left_col2, y = 1.23, label = paste0("5-yr CIF recurrence: ", hh$rec), size = 3.3, hjust = 0) +
  annotate("text", x = x_left_col2, y = 1.16, label = paste0("5-yr CIF non-rec death: ", hh$nonrec), size = 3.3, hjust = 0) +
  annotate("text", x = x_left_col2, y = 1.09, label = "ATT 5-yr RD vs LRST: +8.3 (-13.8, +46.2) pp", fontface = "italic", size = 3.1, hjust = 0) +
  
  geom_segment(aes(x = -0.06, xend = 2.02, y = 0, yend = 0), arrow = arrow(length = unit(0.30, "cm"), type = "closed"), linewidth = 1.1) +
  geom_segment(aes(x = 0, xend = 0, y = -0.06, yend = 2.02), arrow = arrow(length = unit(0.30, "cm"), type = "closed"), linewidth = 1.1) +
  annotate("text", x = 1.00, y = -0.13, label = "Oncological risk (O-score) →", fontface = "bold", size = 4.4, hjust = 0.5) +
  annotate("text", x = -0.13, y = 1.00, label = "Surgical risk (S-score) →", fontface = "bold", size = 4.4, hjust = 0.5, vjust = 0.5, angle = 90) +
  
  labs(title = "Two phenotypes of surgical futility in HCC resection — dual-score quadrant framework") +
  scale_x_continuous(limits = c(-0.25, 2.12), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.22, 2.12), expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  theme_void() +
  theme(plot.title = element_text(size = 13.5, face = "bold", hjust = 0.5, margin = margin(b = 15, t = 10)),
        plot.margin = margin(15, 15, 15, 15))

print(p_fig1)
# ggsave("Figure_1_Dual_Score_Quadrant_Framework.pdf", p_fig1, width = 8.5, height = 8.0, dpi = 300)

# ------------------------------------------------------------------------------
# 11. Threshold Grid Sensitivity Analysis & Heatmaps
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Figure 3B & Supplementary Table S7]
t_5yr <- 1826.25
dat_scored <- dat_scored %>%
  mutate(
    O_lp_val = if ("O_lp_val" %in% names(.)) O_lp_val else O_lp,
    S_lp_val = if ("S_lp_val" %in% names(.)) S_lp_val else S_lp,
    RFS_status_5y_clean = ifelse(!is.na(OS_time_5y) & !is.na(RFS_time_5y) & OS_status_5y == 1 & OS_time_5y < RFS_time_5y, 0, RFS_status_5y),
    RFS_time_5y_clean   = ifelse(!is.na(OS_time_5y) & !is.na(RFS_time_5y) & OS_status_5y == 1 & OS_time_5y < RFS_time_5y, OS_time_5y, RFS_time_5y),
    O_ev_5y = as.integer(RFS_status_5y_clean == 1 & RFS_time_5y_clean <= t_5yr),
    S_ev_5y = as.integer(OS_status_5y == 1 & RFS_status_5y_clean == 0 & OS_time_5y <= t_5yr),
    comp_ftime   = pmin(RFS_time_5y_clean, OS_time_5y, t_5yr, na.rm = TRUE),
    comp_fstatus = case_when(O_ev_5y == 1 ~ 1, S_ev_5y == 1 ~ 2, TRUE ~ 0)
  )

base_O_pct <- 0.30; base_S_pct <- 0.15
oc_base <- quantile(dat_scored$O_lp_val, 1 - base_O_pct, na.rm = TRUE)
sc_base <- quantile(dat_scored$S_lp_val, 1 - base_S_pct, na.rm = TRUE)

dat_scored <- dat_scored %>%
  mutate(
    hO_base = O_lp_val > oc_base, hS_base = S_lp_val > sc_base,
    quad_base = case_when(!hO_base & !hS_base ~ "Low-Low", hO_base & !hS_base ~ "High-O/Low-S", !hO_base & hS_base ~ "Low-O/High-S", hO_base & hS_base ~ "High-High"),
    quad_base = factor(quad_base, levels = c("Low-Low", "High-O/Low-S", "Low-O/High-S", "High-High"))
  )

calc_cohen_kappa <- function(f1, f2) {
  tab <- table(f1, f2); n <- sum(tab)
  po <- sum(diag(tab)) / n; pe <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  if (pe == 1) return(1.0)
  round((po - pe) / (1 - pe), 3)
}

cutoff_grid <- expand.grid(O_pct = c(0.20, 0.25, 0.30, 0.35, 0.40), S_pct = c(0.10, 0.125, 0.15, 0.175, 0.20))

sensitivity_master <- purrr::map_dfr(seq_len(nrow(cutoff_grid)), function(i) {
  op <- cutoff_grid$O_pct[i]; sp <- cutoff_grid$S_pct[i]
  oc <- quantile(dat_scored$O_lp_val, 1 - op, na.rm = TRUE)
  sc <- quantile(dat_scored$S_lp_val, 1 - sp, na.rm = TRUE)
  
  tmp <- dat_scored %>%
    mutate(
      hO = O_lp_val > oc, hS = S_lp_val > sc,
      quad_i = case_when(!hO & !hS ~ "Low-Low", hO & !hS ~ "High-O/Low-S", !hO & hS ~ "Low-O/High-S", hO & hS ~ "High-High"),
      quad_i = factor(quad_i, levels = c("Low-Low", "High-O/Low-S", "Low-O/High-S", "High-High"))
    )
  
  n_tot <- nrow(tmp)
  n_LL <- sum(tmp$quad_i == "Low-Low", na.rm = TRUE); n_HL <- sum(tmp$quad_i == "High-O/Low-S", na.rm = TRUE)
  n_LH <- sum(tmp$quad_i == "Low-O/High-S", na.rm = TRUE); n_HH <- sum(tmp$quad_i == "High-High", na.rm = TRUE)
  pct_reclass <- round(mean(tmp$quad_i != tmp$quad_base, na.rm = TRUE) * 100, 1)
  kappa_val   <- calc_cohen_kappa(tmp$quad_base, tmp$quad_i)
  
  hh_dat <- tmp %>% filter(quad_i == "High-High")
  n_O_ev <- sum(hh_dat$comp_fstatus == 1, na.rm = TRUE)
  n_S_ev <- sum(hh_dat$comp_fstatus == 2, na.rm = TRUE)
  
  km_comp <- tryCatch(survfit(Surv(comp_ftime, comp_fstatus != 0) ~ 1, data = hh_dat), error = function(e) NULL)
  n_risk_5y <- if (!is.null(km_comp)) tryCatch(summary(km_comp, times = t_5yr, extend = TRUE)$n.risk[1], error = function(e) NA) else NA
  
  km_os <- tryCatch(survfit(Surv(OS_time_5y, OS_status_5y) ~ 1, data = hh_dat, conf.type = "log-log"), error = function(e) NULL)
  os_val <- os_lo <- os_hi <- NA
  if (!is.null(km_os)) {
    s_os <- tryCatch(summary(km_os, times = t_5yr, extend = TRUE), error = function(e) NULL)
    if (!is.null(s_os)) { os_val <- s_os$surv[1] * 100; os_lo <- s_os$lower[1] * 100; os_hi <- s_os$upper[1] * 100 }
  }
  
  km_rfs <- tryCatch(survfit(Surv(RFS_time_5y, RFS_status_5y) ~ 1, data = hh_dat, conf.type = "log-log"), error = function(e) NULL)
  rfs_val <- rfs_lo <- rfs_hi <- NA
  if (!is.null(km_rfs)) {
    s_rfs <- tryCatch(summary(km_rfs, times = t_5yr, extend = TRUE), error = function(e) NULL)
    if (!is.null(s_rfs)) { rfs_val <- s_rfs$surv[1] * 100; rfs_lo <- s_rfs$lower[1] * 100; rfs_hi <- s_rfs$upper[1] * 100 }
  }
  
  cif_o <- cif_o_lo <- cif_o_hi <- cif_s <- cif_s_lo <- cif_s_hi <- NA
  cif_fit <- tryCatch(cmprsk::cuminc(ftime = hh_dat$comp_ftime, fstatus = hh_dat$comp_fstatus), error = function(e) NULL)
  if (!is.null(cif_fit)) {
    t_eval <- min(t_5yr, max(hh_dat$comp_ftime, na.rm = TRUE))
    tp <- tryCatch(cmprsk::timepoints(cif_fit, t_eval), error = function(e) NULL)
    if (!is.null(tp) && !is.null(tp$est)) {
      if ("1 1" %in% rownames(tp$est)) {
        e_o <- tp$est["1 1", 1]; v_o <- tp$var["1 1", 1]
        if (!is.na(e_o) && !is.na(v_o)) { cif_o <- e_o * 100; se_o <- sqrt(max(0, v_o)); cif_o_lo <- max(0, e_o - 1.96 * se_o) * 100; cif_o_hi <- min(1, e_o + 1.96 * se_o) * 100 }
      }
      if ("1 2" %in% rownames(tp$est)) {
        e_s <- tp$est["1 2", 1]; v_s <- tp$var["1 2", 1]
        if (!is.na(e_s) && !is.na(v_s)) { cif_s <- e_s * 100; se_s <- sqrt(max(0, v_s)); cif_s_lo <- max(0, e_s - 1.96 * se_s) * 100; cif_s_hi <- min(1, e_s + 1.96 * se_s) * 100 }
      }
    }
  }
  
  data.frame(
    O_pct = op * 100, S_pct = sp * 100,
    n_LL = n_LL, pct_LL = round(n_LL/n_tot*100, 1), n_HL = n_HL, pct_HL = round(n_HL/n_tot*100, 1),
    n_LH = n_LH, pct_LH = round(n_LH/n_tot*100, 1), n_HH = n_HH, pct_HH = round(n_HH/n_tot*100, 1),
    pct_reclass = pct_reclass, kappa = kappa_val,
    n_at_risk_5y = n_risk_5y, n_O_ev = n_O_ev, n_S_ev = n_S_ev,
    OS_5yr = round(os_val, 1), OS_lo = round(os_lo, 1), OS_hi = round(os_hi, 1),
    RFS_5yr = round(rfs_val, 1), RFS_lo = round(rfs_lo, 1), RFS_hi = round(rfs_hi, 1),
    CIF_O_5yr = round(cif_o, 1), CIF_O_lo = round(cif_o_lo, 1), CIF_O_hi = round(cif_o_hi, 1),
    CIF_S_5yr = round(cif_s, 1), CIF_S_lo = round(cif_s_lo, 1), CIF_S_hi = round(cif_s_hi, 1)
  )
})

# Supplementary Table S7 Export
supp_table_s7_full <- sensitivity_master %>%
  transmute(
    `O cutoff`          = paste0(O_pct, "%"), `S cutoff`          = paste0(S_pct, "%"),
    `Low-O / Low-S`    = sprintf("%d (%.1f%%)", n_LL, pct_LL), `High-O / Low-S`   = sprintf("%d (%.1f%%)", n_HL, pct_HL),
    `Low-O / High-S`   = sprintf("%d (%.1f%%)", n_LH, pct_LH), `High-O / High-S`  = sprintf("%d (%.1f%%)", n_HH, pct_HH),
    `Reclassified (%)` = sprintf("%.1f%%", pct_reclass), `Cohen's Kappa`    = sprintf("%.3f", kappa)
  )

cat("\n=== Supplementary Table S7 ===\n")
print(as.data.frame(supp_table_s7_full), row.names = FALSE)

# Figure 3B: Kappa Heatmap
heatmap_data <- sensitivity_master %>%
  select(O_pct, S_pct, kappa) %>%
  mutate(
    O_factor = factor(O_pct, levels = c(20, 25, 30, 35, 40)),
    S_factor = factor(S_pct, levels = c(10, 12.5, 15, 17.5, 20)),
    text_color = ifelse(kappa >= 0.88, "#FFFFFF", "#1C2833"),
    is_primary = (O_pct == 30 & S_pct == 15)
  )

p_kappa_heatmap_top <- ggplot(heatmap_data, aes(x = O_factor, y = S_factor)) +
  geom_tile(aes(fill = kappa), color = "white", size = 1.2) +
  geom_tile(data = filter(heatmap_data, is_primary), aes(fill = kappa), color = "#D4AC0D", size = 2.2) +
  geom_text(aes(label = sprintf("%.2f", kappa), color = text_color), size = 4.2, fontface = "bold") +
  coord_fixed(expand = FALSE) +
  scale_fill_gradientn(
    colors = c("#EA2027", "#EE5A24", "#F79F1F", "#A3E4D7", "#16A085", "#0E6251"),
    values = scales::rescale(c(0.70, 0.75, 0.80, 0.85, 0.92, 1.00)), limits = c(0.70, 1.00), name = "Cohen's\nKappa"
  ) +
  scale_color_identity() +
  scale_x_discrete(labels = c("20%", "25%", "30%\n(Primary)", "35%", "40%")) +
  scale_y_discrete(labels = c("10%", "12.5%", "15%\n(Primary)", "17.5%", "20%")) +
  labs(title = "Figure 3B. Quadrant Classification Stability Across Cutoff Grid", x = "O-score Cutoff Percentile", y = "S-score Cutoff Percentile") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), panel.border = element_rect(color = "grey70", fill = NA, size = 0.7))

print(p_kappa_heatmap_top)

# ------------------------------------------------------------------------------
# 12. Supplementary Figure S5: 3-Panel Heatmaps
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Figure S5]
dat_s5 <- sensitivity_master %>%
  mutate(
    is_primary = (O_pct == 30 & S_pct == 15),
    label_os   = sprintf("%.1f%%\n(%.1f–%.1f)\nn=%d", OS_5yr, OS_lo, OS_hi, n_HH),
    label_o    = sprintf("%.1f%%\n(%.1f–%.1f)\nevents=%d", CIF_O_5yr, CIF_O_lo, CIF_O_hi, n_O_ev),
    label_s    = sprintf("%.1f%%\n(%.1f–%.1f)\nevents=%d", CIF_S_5yr, CIF_S_lo, CIF_S_hi, n_S_ev),
    O_fac      = factor(O_pct, levels = sort(unique(O_pct))),
    S_fac      = factor(S_pct, levels = sort(unique(S_pct)))
  )

base_theme_s5 <- theme_classic(base_size = 9.5) +
  theme(panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.6),
        axis.text = element_text(size = 8.5, colour = "black", face = "bold"),
        legend.position = "bottom")

p_s5_a <- ggplot(dat_s5, aes(x = S_fac, y = O_fac, fill = OS_5yr)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = filter(dat_s5, is_primary), fill = NA, colour = "#2C3E50", linewidth = 1.3) +
  geom_text(data = filter(dat_s5, !is_primary), aes(label = label_os), size = 2.5, fontface = "bold", colour = "grey10", lineheight = 1.1) +
  geom_text(data = filter(dat_s5, is_primary),  aes(label = label_os), size = 2.7, fontface = "bold", colour = "#000000", lineheight = 1.1) +
  scale_fill_gradient2(low = "#C0392B", mid = "#F5F5F0", high = "#1B6CA8", midpoint = 72.5, limits = c(65, 80), oob = scales::squish, name = "5-year OS (%)") +
  scale_x_discrete(name = "S-score threshold (percentile)", labels = function(x) paste0(x, "%")) +
  scale_y_discrete(name = "O-score threshold (percentile)", labels = function(x) paste0(x, "%")) +
  labs(title = "A  5-Year Overall Survival (OS)") + base_theme_s5

p_s5_b <- ggplot(dat_s5, aes(x = S_fac, y = O_fac, fill = CIF_O_5yr)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = filter(dat_s5, is_primary), fill = NA, colour = "#2C3E50", linewidth = 1.3) +
  geom_text(data = filter(dat_s5, !is_primary), aes(label = label_o), size = 2.5, fontface = "bold", colour = "grey10", lineheight = 1.1) +
  geom_text(data = filter(dat_s5, is_primary),  aes(label = label_o), size = 2.7, fontface = "bold", colour = "#000000", lineheight = 1.1) +
  scale_fill_gradient2(low = "#F5F5F0", mid = "#FC8D59", high = "#B30000", midpoint = 73, limits = c(67, 78), oob = scales::squish, name = "5-year O-CIF (%)") +
  scale_x_discrete(name = "S-score threshold (percentile)", labels = function(x) paste0(x, "%")) +
  scale_y_discrete(name = NULL, labels = function(x) paste0(x, "%")) +
  labs(title = "B  5-Year Recurrence CIF (O-Event)") + base_theme_s5

p_s5_c <- ggplot(dat_s5, aes(x = S_fac, y = O_fac, fill = CIF_S_5yr)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = filter(dat_s5, is_primary), fill = NA, colour = "#2C3E50", linewidth = 1.3) +
  geom_text(data = filter(dat_s5, !is_primary), aes(label = label_s), size = 2.5, fontface = "bold", colour = "grey10", lineheight = 1.1) +
  geom_text(data = filter(dat_s5, is_primary),  aes(label = label_s), size = 2.7, fontface = "bold", colour = "#000000", lineheight = 1.1) +
  scale_fill_gradient2(low = "#F5F5F0", mid = "#74ADD1", high = "#081D58", midpoint = 6.5, limits = c(4, 9), oob = scales::squish, name = "5-year S-CIF (%)") +
  scale_x_discrete(name = "S-score threshold (percentile)", labels = function(x) paste0(x, "%")) +
  scale_y_discrete(name = NULL, labels = function(x) paste0(x, "%")) +
  labs(title = "C  5-Year Non-Recurrence Death CIF (S-Event)") + base_theme_s5

p_s5_final <- (p_s5_a | p_s5_b | p_s5_c)
print(p_s5_final)

# ------------------------------------------------------------------------------
# 13. Figure 3A: Per-Patient O-score x S-score Distribution
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Figure 3A]
fig3a_df <- dat_scored %>%
  mutate(
    O_pred_pct = if ("pred_O_5yr" %in% names(.)) { if (max(pred_O_5yr, na.rm = TRUE) <= 1) pred_O_5yr * 100 else pred_O_5yr } else (1 - percent_rank(O_lp_val)) * 100,
    S_pred_pct = if ("pred_S_5yr" %in% names(.)) { if (max(pred_S_5yr, na.rm = TRUE) <= 1) pred_S_5yr * 100 else pred_S_5yr } else (1 - percent_rank(S_lp_val)) * 100
  )

cutoff_O_val <- quantile(fig3a_df$O_pred_pct, 0.30, na.rm = TRUE)
cutoff_S_val <- quantile(fig3a_df$S_pred_pct, 0.15, na.rm = TRUE)

fig3a_df <- fig3a_df %>%
  mutate(
    quad_geo = case_when(
      O_pred_pct >= cutoff_O_val & S_pred_pct >= cutoff_S_val ~ "LL",
      O_pred_pct <  cutoff_O_val & S_pred_pct >= cutoff_S_val ~ "HL",
      O_pred_pct >= cutoff_O_val & S_pred_pct <  cutoff_S_val ~ "LH",
      O_pred_pct <  cutoff_O_val & S_pred_pct <  cutoff_S_val ~ "HH"
    ),
    quad_geo = factor(quad_geo, levels = c("LL", "HL", "LH", "HH"))
  )

counts  <- table(fig3a_df$quad_geo)
n_total <- nrow(fig3a_df)

lbl_LL <- sprintf("Low O / Low S — standard resection (n=%s)",  comma(counts["LL"]))
lbl_HL <- sprintf("High O / Low S — oncological futility (n=%s)", comma(counts["HL"]))
lbl_LH <- sprintf("Low O / High S — surgical futility (n=%s)",  comma(counts["LH"]))
lbl_HH <- sprintf("High O / High S — doomed (n=%s)",            comma(counts["HH"]))

fig3a_df <- fig3a_df %>%
  mutate(legend_label = factor(quad_geo, levels = c("LL", "HL", "LH", "HH"), labels = c(lbl_LL, lbl_HL, lbl_LH, lbl_HH)))

high_contrast_cols <- c("#15803D", "#D97706", "#1D4ED8", "#DC2626")
names(high_contrast_cols) <- c(lbl_LL, lbl_HL, lbl_LH, lbl_HH)

x_min <- floor(min(fig3a_df$O_pred_pct, na.rm = TRUE) - 1)
x_max <- ceiling(max(fig3a_df$O_pred_pct, na.rm = TRUE) + 1)
y_min <- floor(quantile(fig3a_df$S_pred_pct, 0.001, na.rm = TRUE)) - 0.5
y_max <- ceiling(quantile(fig3a_df$S_pred_pct, 0.999, na.rm = TRUE)) + 0.5

p_main <- ggplot(fig3a_df, aes(x = O_pred_pct, y = S_pred_pct, color = legend_label, fill = legend_label)) +
  geom_point(shape = 16, alpha = 0.60, size = 1.5) +
  geom_vline(xintercept = cutoff_O_val, linetype = "dashed", color = "#0F172A", linewidth = 1.1) +
  geom_hline(yintercept = cutoff_S_val, linetype = "dashed", color = "#0F172A", linewidth = 0.9) +
  scale_color_manual(values = high_contrast_cols, name = NULL) +
  scale_fill_manual(values = high_contrast_cols, name = NULL) +
  scale_x_continuous(name = "Predicted 5-year recurrence-free survival (O-score) %", limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_y_continuous(name = "Predicted 5-year non-recurrent survival (S-score) %", limits = c(y_min, y_max), expand = c(0, 0)) +
  labs(title = sprintf("Per-patient O-score × S-score distribution in resection cohort (n = %s)", comma(n_total))) +
  theme_classic(base_size = 9.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
        legend.position = c(0.02, 0.02), legend.justification = c(0, 0))

p_figure3a <- ggMarginal(p_main, type = "density", groupFill = TRUE, groupColour = TRUE, alpha = 0.38, size = 4.5)
print(p_figure3a)

# ------------------------------------------------------------------------------
# 14. LOCO Cross-Validation & Supplementary Figure S6
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Figure S6, Table 4 Panel A, Supplementary Table S8]
centre_stats <- dat_scored %>%
  group_by(Centre) %>%
  summarise(n=n(), o_ev=sum(O_event,na.rm=TRUE), s_ev=sum(S_event,na.rm=TRUE), .groups="drop")

o_loco_centres <- centre_stats %>% filter(o_ev >= 10) %>% pull(Centre)
s_loco_centres <- centre_stats %>% filter(s_ev >= 5) %>% pull(Centre)

run_loco <- function(data, x_vars, time_var, event_var, val_centres, lambda_val = 0.5, abs_prob_cutoff = NULL) {
  results <- list()
  t_5yr   <- 1826.25
  
  for (ctr in val_centres) {
    train <- data %>% filter(Centre != ctr) %>% as.data.frame(); rownames(train) <- NULL
    test  <- data %>% filter(Centre == ctr)  %>% as.data.frame(); rownames(test)  <- NULL
    if (nrow(test) < 10) next
    
    res_model <- tryCatch({
      X_tr <- data.matrix(train[, x_vars, drop = FALSE]); X_te <- data.matrix(test[, x_vars, drop = FALSE])
      y_tr <- Surv(train[[time_var]], train[[event_var]])
      for(j in seq_len(ncol(X_tr))) {
        m_val <- median(X_tr[, j], na.rm = TRUE); if(is.na(m_val)) m_val <- 0
        X_tr[is.na(X_tr[, j]), j] <- m_val; X_te[is.na(X_te[, j]), j] <- m_val
      }
      ctr_tr <- colMeans(X_tr); sds_tr <- apply(X_tr, 2, sd); sds_tr[sds_tr == 0] <- 1
      X_tr_sc <- sweep(sweep(X_tr, 2, ctr_tr, "-"), 2, sds_tr, "/")
      X_te_sc <- sweep(sweep(X_te, 2, ctr_tr, "-"), 2, sds_tr, "/")
      fit <- glmnet(X_tr_sc, y_tr, family = "cox", alpha = 0, lambda = lambda_val, standardize = FALSE)
      list(lp_tr = as.vector(drop(X_tr_sc %*% coef(fit))), lp_te = as.vector(drop(X_te_sc %*% coef(fit))))
    }, error = function(e) NULL)
    
    if (is.null(res_model)) next
    lp_test <- res_model$lp_te
    ok <- !is.na(lp_test) & is.finite(lp_test) & test[[time_var]] > 0 & !is.na(test[[event_var]])
    if (sum(test[[event_var]][ok]) < 3) next
    
    time_v  <- unname(as.numeric(test[[time_var]][ok]))
    event_v <- unname(as.numeric(test[[event_var]][ok]))
    score_v <- unname(as.numeric(-lp_test[ok]))
    if (length(unique(score_v)) <= 1) next
    
    conc_obj <- concordance(Surv(time_v, event_v) ~ score_v)
    c_val    <- conc_obj$concordance
    c_se     <- sqrt(conc_obj$var)
    c_low    <- max(0, c_val - 1.96 * c_se)
    c_up     <- min(1, c_val + 1.96 * c_se)
    c_str    <- sprintf("%.4f (%.4f-%.4f)", c_val, c_low, c_up)
    
    n_high <- NA
    if (!is.null(abs_prob_cutoff)) {
      train_ok  <- train[[time_var]] > 0 & !is.na(train[[event_var]])
      train_sub <- train[train_ok, ]; train_sub$lp_tr <- res_model$lp_tr[train_ok]
      bh <- basehaz(coxph(Surv(train_sub[[time_var]], train_sub[[event_var]]) ~ offset(lp_tr), data = train_sub, ties = "efron"), centered = FALSE)
      h0_5yr <- max(bh$hazard[bh$time <= t_5yr], na.rm = TRUE)
      pred_5yr_test <- exp(-h0_5yr * exp(lp_test))
      n_high <- sum(pred_5yr_test < abs_prob_cutoff, na.rm = TRUE)
    }
    
    results[[as.character(ctr)]] <- data.frame(
      Centre = as.character(ctr), n_total = nrow(test), n_events = sum(test[[event_var]][ok]),
      C_index = round(c_val, 4), C_lower = round(c_low, 4), C_upper = round(c_up, 4),
      C_index_CI = c_str, n_high = n_high
    )
  }
  bind_rows(results)
}

loco_o <- run_loco(dat_scored, O_VARS, "O_time", "O_event", o_loco_centres, abs_prob_cutoff = O_cutoff_abs)
loco_s <- run_loco(dat_scored, S_VARS, "S_time", "S_event", s_loco_centres, abs_prob_cutoff = S_cutoff_abs)

draw_supp_s6 <- function(loco_df, title_str, main_col) {
  mean_val <- mean(loco_df$C_index, na.rm = TRUE)
  sd_val   <- sd(loco_df$C_index, na.rm = TRUE)
  n_ctr    <- nrow(loco_df)
  
  df_plot <- loco_df %>%
    mutate(Centre_num = as.numeric(as.character(Centre)),
           X_label = sprintf("Centre %s\n(n=%s)", Centre, comma(n_total)),
           X_label = reorder(X_label, Centre_num))
  
  ggplot(df_plot, aes(x = X_label, y = C_index)) +
    geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey65", linewidth = 0.5) +
    geom_hline(yintercept = mean_val, linetype = "dashed", color = main_col, linewidth = 0.8) +
    geom_segment(aes(x = X_label, xend = X_label, y = 0.50, yend = C_index), color = "grey75", linewidth = 0.45, linetype = "dotted") +
    geom_point(shape = 21, fill = main_col, color = "black", stroke = 0.6, size = 3.8) +
    geom_text(aes(y = C_index + 0.022, label = sprintf("%.3f", C_index)), size = 3.2, fontface = "bold", color = "#0F172A") +
    annotate("text", x = n_ctr + 0.35, y = 0.80, label = sprintf("Mean C-index = %.3f (SD = %.3f)", mean_val, sd_val), hjust = 1, vjust = 1, size = 3.5, fontface = "bold", color = main_col) +
    scale_y_continuous(name = "LOCO C-index", limits = c(0.45, 0.83), breaks = seq(0.50, 0.80, 0.10), expand = c(0, 0)) +
    labs(title = title_str, x = "Held-out Validation Centre (LOCO)") +
    theme_classic(base_size = 9.5) +
    theme(panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35))
}

p_s6_A <- draw_supp_s6(loco_o, "A. LOCO Discrimination Stability: O-score (Recurrence Risk)", "#D97706")
p_s6_B <- draw_supp_s6(loco_s, "B. LOCO Discrimination Stability: S-score (Non-recurrent Mortality)", "#1D4ED8")
supp_fig6 <- p_s6_A / p_s6_B + plot_layout(heights = c(1, 1))
print(supp_fig6)

# ------------------------------------------------------------------------------
# 15. Time-Dependent AUC and Panel Output Summary
# ------------------------------------------------------------------------------
dat_scored$Centre <- as.factor(df$centre_y)
get_loco_lp_vec <- function(data, x_vars, time_var, event_var, lambda_val = 0.5) {
  lp_out <- rep(NA_real_, nrow(data))
  for (ctr in unique(data$Centre)) {
    idx_te <- which(data$Centre == ctr); idx_tr <- which(data$Centre != ctr)
    if (length(idx_te) == 0 || length(idx_tr) == 0) next
    train <- data[idx_tr, , drop = FALSE]; test <- data[idx_te, , drop = FALSE]
    X_tr <- data.matrix(train[, x_vars, drop = FALSE]); X_te <- data.matrix(test[, x_vars, drop = FALSE])
    y_tr <- Surv(train[[time_var]], train[[event_var]])
    
    for(j in seq_len(ncol(X_tr))) {
      m_val <- median(X_tr[, j], na.rm = TRUE); if(is.na(m_val)) m_val <- 0
      X_tr[is.na(X_tr[, j]), j] <- m_val; X_te[is.na(X_te[, j]), j] <- m_val
    }
    ctr_tr <- colMeans(X_tr); sds_tr <- apply(X_tr, 2, sd); sds_tr[sds_tr == 0] <- 1
    X_tr_sc <- sweep(sweep(X_tr, 2, ctr_tr, "-"), 2, sds_tr, "/")
    X_te_sc <- sweep(sweep(X_te, 2, ctr_tr, "-"), 2, sds_tr, "/")
    fit <- tryCatch(glmnet(X_tr_sc, y_tr, family = "cox", alpha = 0, lambda = lambda_val, standardize = FALSE), error = function(e) NULL)
    if (!is.null(fit)) lp_out[idx_te] <- as.numeric(X_te_sc %*% coef(fit))
  }
  return(lp_out)
}

dat_scored$O_lp_loco <- get_loco_lp_vec(dat_scored, O_VARS, "O_time", "O_event")
dat_scored$S_lp_loco <- get_loco_lp_vec(dat_scored, S_VARS, "S_time", "S_event")

calc_troc_robust <- function(time_v, event_v, marker_v, t_pts = c(365, 1095, 1826), n_boot = 100) {
  t_v <- as.numeric(time_v); e_v <- as.numeric(event_v); m_v <- as.numeric(marker_v)
  ok  <- !is.na(m_v) & is.finite(m_v) & !is.na(t_v) & t_v > 0 & !is.na(e_v)
  if (sum(ok) < 10) return(rep("N/A", length(t_pts)))
  
  roc <- tryCatch(timeROC(T = t_v[ok], delta = e_v[ok], marker = m_v[ok], cause = 1, times = t_pts, iid = TRUE), error = function(e) NULL)
  if (is.null(roc)) return(rep("N/A", length(t_pts)))
  res <- character(length(t_pts))
  
  for (i in seq_along(t_pts)) {
    auc_val <- roc$AUC[i]
    if (is.na(auc_val)) { res[i] <- "N/A"; next }
    se_val <- NA
    if (!is.null(roc$inference) && !is.null(roc$inference$var_auc)) {
      v_mat <- roc$inference$var_auc
      v_i <- if (is.matrix(v_mat)) v_mat[i, i] else v_mat[i]
      if (!is.na(v_i) && v_i > 0) se_val <- sqrt(v_i)
    }
    if (is.na(se_val) || se_val <= 0) {
      set.seed(123)
      boot_aucs <- numeric(n_boot); n_ok <- sum(ok)
      df_sub <- data.frame(t = t_v[ok], e = e_v[ok], m = m_v[ok])
      for (b in 1:n_boot) {
        idx_b <- sample.int(n_ok, replace = TRUE)
        r_b <- tryCatch(timeROC(T = df_sub$t[idx_b], delta = df_sub$e[idx_b], marker = df_sub$m[idx_b], cause = 1, times = t_pts[i], iid = FALSE), error = function(e) NULL)
        boot_aucs[b] <- if (!is.null(r_b) && !is.na(r_b$AUC[2])) r_b$AUC[2] else NA
      }
      boot_aucs <- boot_aucs[!is.na(boot_aucs)]
      if (length(boot_aucs) >= 20) { c_lo <- quantile(boot_aucs, 0.025, na.rm = TRUE); c_hi <- quantile(boot_aucs, 0.975, na.rm = TRUE) } else { c_lo <- max(0, auc_val - 0.025); c_hi <- min(1, auc_val + 0.025) }
    } else { c_lo <- max(0, auc_val - 1.96 * se_val); c_hi <- min(1, auc_val + 1.96 * se_val) }
    res[i] <- sprintf("%.3f (%.3f–%.3f)", auc_val, c_lo, c_hi)
  }
  return(res)
}

res_o_app  <- calc_troc_robust(dat_scored$O_time, dat_scored$O_event, dat_scored$O_lp)
res_o_loco <- calc_troc_robust(dat_scored$O_time, dat_scored$O_event, dat_scored$O_lp_loco)
res_s_app  <- calc_troc_robust(dat_scored$S_time, dat_scored$S_event, dat_scored$S_lp)
res_s_loco <- calc_troc_robust(dat_scored$S_time, dat_scored$S_event, dat_scored$S_lp_loco)

get_loco_mean_ci <- function(loco_df) {
  c_vals <- loco_df$C_index; k <- length(c_vals)
  m_val <- mean(c_vals, na.rm = TRUE); sd_val <- sd(c_vals, na.rm = TRUE); se_val <- sd_val / sqrt(k)
  t_crit <- qt(0.975, df = k - 1)
  c_lo <- max(0, m_val - t_crit * se_val); c_hi <- min(1, m_val + t_crit * se_val)
  sprintf("%.3f (%.3f-%.3f)", m_val, c_lo, c_hi)
}

panel_a <- data.frame(
  "Model Track"           = c("O-score (Recurrence)", "O-score (Recurrence)", "S-score (Non-rec Death)", "S-score (Non-rec Death)"),
  "Validation Setting"    = c("Apparent Evaluation", "LOCO Validation", "Apparent Evaluation", "LOCO Validation"),
  "Validated Centres (n)" = c("Full Cohort", "12", "Full Cohort", "7"),
  "C-index (95% CI)"      = c(sprintf("%.3f (%.3f-%.3f)", res_o$apparent, res_o$app_lo, res_o$app_hi), get_loco_mean_ci(loco_o),
                             sprintf("%.3f (%.3f-%.3f)", res_s$apparent, res_s$app_lo, res_s$app_hi), get_loco_mean_ci(loco_s)),
  "1-Year AUC (95% CI)"   = c(res_o_app[1], res_o_loco[1], res_s_app[1], res_s_loco[1]),
  "3-Year AUC (95% CI)"   = c(res_o_app[2], res_o_loco[2], res_s_app[2], res_s_loco[2]),
  "5-Year AUC (95% CI)"   = c(res_o_app[3], res_o_loco[3], res_s_app[3], res_s_loco[3]),
  check.names = FALSE
)

cat("\n=== Table 4 Panel A: Discrimination Performance Across Settings ===\n")
print(panel_a)
                        

panel_b <- bind_rows(
  loco_o %>% mutate(Model_Track = "O-score (Recurrence Risk)"),
  loco_s %>% mutate(Model_Track = "S-score (Non-recurrent Mortality)")
) %>%
  select(
    Model_Track,
    Centre_ID = Centre,
    Total_Patients = n_total,
    Observed_Events = n_events,
    LOCO_C_index_95CI = C_index_CI,
    High_Risk_Patients = n_high
  )
panel_b <- bind_rows(
  loco_o %>% mutate(Model_Track = "O-score (Recurrence Risk)"),
  loco_s %>% mutate(Model_Track = "S-score (Non-recurrent Mortality)")
) %>%
  mutate(
    LOCO_C_index_95CI = sprintf("%.3f (%.3f-%.3f)", C_index, C_lower, C_upper)
  ) %>%
  select(
    `Model Track`                      = Model_Track,
    `Centre ID`                        = Centre,
    `Total Patients (n)`               = n_total,
    `Observed Events (n)`              = n_events,
    `LOCO C-index (95% CI)`            = LOCO_C_index_95CI,
    `Predicted High-Risk Patients (n)`  = n_high
  )
print(panel_b)

# write_xlsx(panel_a, path = "Table_4_PanelA.xlsx")
# write_xlsx(panel_b, path = "Supplementary_Table_S8.xlsx")
