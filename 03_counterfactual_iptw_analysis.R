# ==============================================================================
# Script Name: 03_counterfactual_iptw_analysis.R
# Purpose: LRST comparator imputation, O/S risk score projection, IPTW-ATT weighting, 
#          covariate balance diagnostics, doubly robust multivariable Cox models, 
#          and sensitivity analyses across risk quadrants.
# Outputs: Table 3, Figure 2, Supplementary Figures S3, S4,
#          Supplementary Tables S3, S4, S5, S6.
# Language: R (>= 4.0.0)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(survival)
library(glmnet)
library(mice)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(purrr)
library(tidyr)

# ------------------------------------------------------------------------------
# 2. Ensure Resection Data State (dat_resect)
# ------------------------------------------------------------------------------
# Retain pristine resection cohort to prevent workspace pollution
dat_resect <- dat_scored

if (exists("df")) {
  add_cols <- intersect(c("O_time", "O_event", "S_time", "S_event", "age", "gender", "DM", "HTN", "nash", "centre_x"), colnames(df))
  for (col in add_cols) {
    if (!col %in% colnames(dat_resect)) dat_resect[[col]] <- df[[col]]
  }
}

# ------------------------------------------------------------------------------
# 3. Import and Preprocess LRST Comparator Cohort
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
  mutate(across(everything(), ~ ifelse(. %in% c(".", "NaN", " ", ""), NA, .))) %>%
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

# Independent MICE Imputation for LRST Cohort (m = 20)
vars_imp <- unique(c(O_VARS, S_VARS, "OS_time_5y", "OS_status_5y", "sn", "centre_x", "TACE", "Chemo", "RE", "age", "gender", "DM", "HTN", "nash"))
df_imp_in <- df_LRST %>% select(all_of(intersect(vars_imp, colnames(df_LRST))))

imp_method_lrst <- make.method(df_imp_in)
imp_method_lrst[intersect(c("OS_time_5y", "OS_status_5y", "sn", "centre_x"), names(imp_method_lrst))] <- ""

cat(">> Running MICE imputation for LRST cohort (m = 20)...\n")
imp_obj_lrst   <- mice(df_imp_in, m = 20, maxit = 10, method = imp_method_lrst, seed = 2026, printFlag = FALSE)
df_LRST_scored <- complete(imp_obj_lrst, 1)

df_LRST_scored$OS_time_5y   <- df_imp_in$OS_time_5y
df_LRST_scored$OS_status_5y <- df_imp_in$OS_status_5y
df_LRST_scored$sn           <- df_imp_in$sn
df_LRST_scored$First_ever   <- df_LRST$First_ever

df_LRST_scored$RFS_time_5y   <- suppressWarnings(as.numeric(df_LRST$RFS_time_5y))
df_LRST_scored$RFS_status_5y <- suppressWarnings(as.numeric(df_LRST$RFS_status_5y))

# ------------------------------------------------------------------------------
# 4. Project Resection-Trained Ridge Cox Scores onto LRST Cohort
# ------------------------------------------------------------------------------
df_LRST_scored$O_lp <- get_lp_ridge(train_dat = dat_resect, test_dat = df_LRST_scored, x_vars = O_VARS, time_var = "O_time", event_var = "O_event", lambda_val = 0.5)$lp
df_LRST_scored$S_lp <- get_lp_ridge(train_dat = dat_resect, test_dat = df_LRST_scored, x_vars = S_VARS, time_var = "S_time", event_var = "S_event", lambda_val = 0.5)$lp

fit_O_lp <- coxph(Surv(O_time, O_event) ~ offset(O_lp), data = dat_resect)
fit_S_lp <- coxph(Surv(S_time, S_event) ~ offset(S_lp), data = dat_resect)

t_5yr <- ifelse(max(dat_resect$OS_time_5y, na.rm = TRUE) > 100, 1826.25, 5)
bh_O  <- basehaz(fit_O_lp, centered = FALSE); H0_O <- bh_O$hazard[max(which(bh_O$time <= t_5yr))]
bh_S  <- basehaz(fit_S_lp, centered = FALSE); H0_S <- bh_S$hazard[max(which(bh_S$time <= t_5yr))]

df_LRST_scored$pred_O_5yr <- exp(-H0_O * exp(df_LRST_scored$O_lp))
df_LRST_scored$pred_S_5yr <- exp(-H0_S * exp(df_LRST_scored$S_lp))

df_LRST_scored <- df_LRST_scored %>%
  mutate(
    high_O = as.integer(pred_O_5yr < O_cutoff_abs),
    high_S = as.integer(pred_S_5yr < S_cutoff_abs),
    quadrant = case_when(
      high_O == 0 & high_S == 0 ~ "Low-O / Low-S",
      high_O == 1 & high_S == 0 ~ "High-O / Low-S",
      high_O == 0 & high_S == 1 ~ "Low-O / High-S",
      high_O == 1 & high_S == 1 ~ "High-O / High-S",
      TRUE ~ "Low-O / Low-S"
    ),
    quadrant = factor(quadrant, levels = c("Low-O / Low-S", "High-O / Low-S", "Low-O / High-S", "High-O / High-S")),
    LRST_type = case_when(
      "TACE" %in% colnames(.) & TACE == 1 ~ "TACE",
      "Chemo" %in% colnames(.) & Chemo == 1 ~ "Systemic",
      "RE" %in% colnames(.) & RE == 1 ~ "RE",
      TRUE ~ "TACE"
    )
  )

# Combine Resection (trt=1) and LRST (trt=0) into dat_pooled
if (!"sn" %in% colnames(dat_resect)) dat_resect$sn <- paste0("RESECT_", 1:nrow(dat_resect))
if (!"LRST_type" %in% colnames(dat_resect)) dat_resect$LRST_type <- "Resection"

common_cols <- intersect(colnames(dat_resect), colnames(df_LRST_scored))
d_resect    <- dat_resect      %>% mutate(trt = 1) %>% select(trt, all_of(common_cols))
d_lrst      <- df_LRST_scored  %>% mutate(trt = 0) %>% select(trt, all_of(common_cols))

dat_pooled  <- bind_rows(d_resect, d_lrst) %>% filter(!is.na(OS_time_5y) & !is.na(OS_status_5y))

cat("\n>> Treatment group distribution:\n")
cat(sprintf("   Resection (trt=1): n = %d\n", sum(dat_pooled$trt == 1, na.rm=TRUE)))
cat(sprintf("   LRST (trt=0):      n = %d\n", sum(dat_pooled$trt == 0, na.rm=TRUE)))
cat("\n>> Quadrant x Treatment cross-tabulation:\n")
print(table(dat_pooled$quadrant, dat_pooled$trt, dnn = c("Quadrant", "Treatment (1=Resection)")))

# ------------------------------------------------------------------------------
# 5. Propensity Score Estimation and ATT Weighting
# ------------------------------------------------------------------------------
PS_VARS <- unique(c(O_VARS, S_VARS, "age", "gender", "nash", "etioletoh"))
ps_vars_in_data <- intersect(PS_VARS, colnames(dat_pooled))

ps_formula <- as.formula(paste("trt ~", paste(ps_vars_in_data, collapse = " + ")))
ps_model <- glm(ps_formula, data = dat_pooled, family = binomial(link = "logit"))
dat_pooled$ps <- predict(ps_model, newdata = dat_pooled, type = "response")

# Calculate Stabilized and Winsorized ATT Weights
dat_pooled <- dat_pooled %>%
  mutate(
    ps_trim   = pmin(pmax(ps, 0.01), 0.99),
    w_ATT     = ifelse(trt == 1, 1.0, ps_trim / (1 - ps_trim)),
    w_ATT_win = pmin(w_ATT, quantile(w_ATT, 0.99, na.rm=TRUE))
  )

# ------------------------------------------------------------------------------
# 6. Supplementary Figure S3: Covariate Balance (Dumbbell Love Plot)
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Figure S3]
smd_calc <- function(var, trt, wt = NULL) {
  x1 <- var[trt==1 & !is.na(var)]; x0 <- var[trt==0 & !is.na(var)]
  if (is.null(wt)) {
    d <- (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2)
  } else {
    w1 <- wt[trt==1 & !is.na(var)]; w0 <- wt[trt==0 & !is.na(var)]
    m1 <- weighted.mean(x1, w1); m0 <- weighted.mean(x0, w0)
    v1 <- sum(w1*(x1-m1)^2)/sum(w1); v0 <- sum(w0*(x0-m0)^2)/sum(w0)
    d <- (m1 - m0) / sqrt((v1 + v0) / 2)
  }
  return(d)
}

balance_tbl <- map_dfr(ps_vars_in_data, function(v) {
  x <- as.numeric(dat_pooled[[v]])
  data.frame(
    Variable   = v,
    SMD_before = round(smd_calc(x, dat_pooled$trt), 3),
    SMD_after  = round(smd_calc(x, dat_pooled$trt, dat_pooled$w_ATT_win), 3)
  )
}) %>% mutate(Balanced = abs(SMD_after) < 0.10)

balance_ordered <- balance_tbl %>%
  mutate(Var_Clean = ifelse(Variable %in% names(var_labels), var_labels[Variable], Variable),
         abs_before = abs(SMD_before)) %>%
  arrange(abs_before) %>%
  mutate(Var_Clean = factor(Var_Clean, levels = Var_Clean))

balance_long <- balance_ordered %>%
  pivot_longer(cols = c(SMD_before, SMD_after), names_to = "Time", values_to = "SMD") %>%
  mutate(Time = ifelse(Time == "SMD_before", "Before IPTW", "After IPTW"),
         Time = factor(Time, levels = c("Before IPTW", "After IPTW")))

p_love_top <- ggplot() +
  annotate("rect", xmin = -0.10, xmax = 0.10, ymin = 0.4, ymax = length(unique(balance_long$Var_Clean)) + 0.6, fill = "#34495E", alpha = 0.08) +
  geom_segment(data = balance_ordered, aes(x = SMD_before, xend = SMD_after, y = Var_Clean, yend = Var_Clean), color = "grey50", size = 0.75) +
  geom_vline(xintercept = 0, color = "black", size = 0.55) +
  geom_vline(xintercept = c(-0.10, 0.10), linetype = "dashed", color = "grey35", size = 0.45) +
  geom_point(data = balance_long, aes(x = SMD, y = Var_Clean, color = Time, shape = Time, fill = Time), size = 3.6, stroke = 0.7) +
  scale_color_manual(values = c("Before IPTW" = "#DF4F49", "After IPTW" = "#00A087")) +
  scale_fill_manual(values = c("Before IPTW" = "#DF4F49", "After IPTW" = "#00A087")) +
  scale_shape_manual(values = c("Before IPTW" = 21, "After IPTW" = 23)) +
  scale_x_continuous(breaks = seq(-0.8, 1.0, by = 0.2), limits = c(-0.85, 1.05), expand = c(0, 0)) +
  labs(title = "Supplementary Figure S3. Assessment of Covariate Balance Before and After IPTW",
       subtitle = "Standardised mean differences (SMDs) for baseline covariates in unweighted and weighted pseudo-populations",
       x = "Standardised Mean Difference (SMD)", y = NULL, color = NULL, shape = NULL, fill = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

print(p_love_top)

# ------------------------------------------------------------------------------
# 7. Stratified Bootstrap for ATT 5-Year Risk Differences (1000 Replicates)
# ------------------------------------------------------------------------------
calc_ess <- function(weights) {
  w <- weights[!is.na(weights) & weights > 0]
  if (length(w) == 0) return(0)
  round((sum(w)^2) / sum(w^2), 1)
}

calc_exact_os_diff <- function(df) {
  quads <- levels(df$quadrant)
  res_list <- list()
  for (q in quads) {
    sub_df <- df %>% filter(quadrant == q & !is.na(w_ATT_win) & !is.na(OS_time_5y))
    get_os_and_ess <- function(trt_val) {
      d_sub <- sub_df %>% filter(trt == trt_val)
      n_raw <- nrow(d_sub)
      if (n_raw < 2) return(list(os = NA_real_, n = n_raw, ess = 0))
      ess_val <- calc_ess(d_sub$w_ATT_win)
      km <- tryCatch(survfit(Surv(OS_time_5y, OS_status_5y) ~ 1, weights = w_ATT_win, data = d_sub), error = function(e) NULL)
      if (is.null(km)) return(list(os = NA_real_, n = n_raw, ess = ess_val))
      s_sum <- summary(km, times = 1826.25, extend = TRUE)
      os_val <- if (length(s_sum$surv) > 0) s_sum$surv[1] * 100 else NA_real_
      list(os = os_val, n = n_raw, ess = ess_val)
    }
    r_res <- get_os_and_ess(1); l_res <- get_os_and_ess(0)
    res_list[[q]] <- data.frame(
      quadrant = q, n_resect = r_res$n, ess_resect = r_res$ess, n_lrst = l_res$n, ess_lrst = l_res$ess,
      ATT_OS_resect = r_res$os, ATT_OS_lrst = l_res$os,
      OS_diff = if (!is.na(r_res$os) && !is.na(l_res$os)) r_res$os - l_res$os else NA_real_
    )
  }
  bind_rows(res_list)
}

os_diff_point <- calc_exact_os_diff(dat_pooled)

cat(">> Running Stratified Bootstrap-1000 for ATT Risk Difference...\n")
B <- 1000
set.seed(2026)
quads <- levels(dat_pooled$quadrant)
boot_mat <- matrix(NA, nrow = B, ncol = length(quads))
colnames(boot_mat) <- quads

for (b in seq_len(B)) {
  boot_df <- dat_pooled %>% group_by(quadrant, trt) %>% sample_frac(replace = TRUE) %>% ungroup()
  ps_b <- tryCatch({ suppressWarnings({ m_b <- glm(ps_formula, data = boot_df, family = binomial(link = "logit")); predict(m_b, newdata = boot_df, type = "response") }) }, error = function(e) boot_df$ps)
  ps_b_clean        <- pmin(pmax(ps_b, 0.01), 0.99)
  w_ATT_b_raw       <- ifelse(boot_df$trt == 1, 1.0, ps_b_clean / (1 - ps_b_clean))
  boot_df$w_ATT_win <- pmin(w_ATT_b_raw, quantile(w_ATT_b_raw, 0.99, na.rm = TRUE))
  res_b <- tryCatch(calc_exact_os_diff(boot_df), error = function(e) NULL)
  if (!is.null(res_b) && nrow(res_b) > 0) boot_mat[b, res_b$quadrant] <- res_b$OS_diff
}

fig1_results <- os_diff_point %>%
  mutate(
    CI_lo = apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)[quadrant],
    CI_hi = apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)[quadrant],
    ATT_OS_resect = round(ATT_OS_resect, 1), ATT_OS_lrst = round(ATT_OS_lrst, 1),
    OS_diff_pp = round(OS_diff, 1), CI_lo = round(CI_lo, 1), CI_hi = round(CI_hi, 1),
    `Figure1_Bottom_Text` = sprintf("ATT 5-yr RD vs LRST: %+.1f (%+.1f, %+.1f) pp", OS_diff_pp, CI_lo, CI_hi)
  )

cat("\n=== Figure 1 Data Output Summary ===\n")
print(fig1_results, row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. Supplementary Table S6: Doubly Adjusted Sensitivity Analysis
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Table S6]
imbalanced_vars <- c("log_Diameter", "HBsAg", "cirr", "logAFP_plus_1")
adj_vars <- intersect(imbalanced_vars, colnames(dat_pooled))

full_adj_formula <- as.formula(paste("Surv(OS_time_5y, OS_status_5y) ~ trt * quadrant +", paste(adj_vars, collapse = " + ")))
fit_full_adj     <- coxph(full_adj_formula, data = dat_pooled, weights = w_ATT_win, robust = TRUE)

reduced_formula  <- as.formula(paste("Surv(OS_time_5y, OS_status_5y) ~ trt + quadrant +", paste(adj_vars, collapse = " + ")))
fit_reduced      <- coxph(reduced_formula, data = dat_pooled, weights = w_ATT_win, robust = TRUE)

lrt_stat    <- 2 * (fit_full_adj$loglik[2] - fit_reduced$loglik[2])
df_diff     <- length(fit_full_adj$coefficients) - length(fit_reduced$coefficients)
p_inter_val <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
p_inter_str <- ifelse(!is.na(p_inter_val) && p_inter_val < 0.001, "< 0.001", sprintf("%.3f", p_inter_val))

s6_table <- data.frame(
  Quadrant = quads, N_Patients = NA_character_, IPTW_RD_95CI = NA_character_,
  Outcome_Adj_RD = NA_character_, Outcome_Adj_HR_95CI = NA_character_, P_value = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(quads)) {
  q_name <- quads[i]
  q_df   <- dat_pooled %>% filter(quadrant == q_name & !is.na(w_ATT_win))
  
  s6_table$N_Patients[i] <- sprintf("%d / %d", sum(q_df$trt == 1, na.rm = TRUE), sum(q_df$trt == 0, na.rm = TRUE))
  
  f1_row <- fig1_results %>% filter(quadrant == q_name)
  if (nrow(f1_row) > 0) { s6_table$IPTW_RD_95CI[i] <- sprintf("%+.1f pp (%+.1f to %+.1f)", f1_row$OS_diff_pp, f1_row$CI_lo, f1_row$CI_hi) }
  
  q_df_resect <- q_df; q_df_resect$trt <- 1
  q_df_lrst   <- q_df; q_df_lrst$trt   <- 0
  pred_adj_resect <- mean(pmax(0, pmin(1, summary(survfit(fit_full_adj, newdata = q_df_resect), times = t_5yr)$surv)))
  pred_adj_lrst   <- mean(pmax(0, pmin(1, summary(survfit(fit_full_adj, newdata = q_df_lrst), times = t_5yr)$surv)))
  s6_table$Outcome_Adj_RD[i] <- sprintf("%+.1f pp", (pred_adj_resect - pred_adj_lrst) * 100)
  
  adj_quad_formula <- as.formula(paste("Surv(OS_time_5y, OS_status_5y) ~ trt +", paste(adj_vars, collapse = " + ")))
  fit_q_adj        <- coxph(adj_quad_formula, data = q_df, weights = w_ATT_win, robust = TRUE)
  sum_q_adj        <- summary(fit_q_adj)
  
  hr_val <- sum_q_adj$coefficients["trt", "exp(coef)"]
  hr_lo  <- sum_q_adj$conf.int["trt", "lower .95"]
  hr_hi  <- sum_q_adj$conf.int["trt", "upper .95"]
  p_val  <- sum_q_adj$coefficients["trt", "Pr(>|z|)"]
  
  s6_table$Outcome_Adj_HR_95CI[i] <- sprintf("%.2f (%.2f–%.2f)", hr_val, hr_lo, hr_hi)
  s6_table$P_value[i]             <- ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val))
}

colnames(s6_table) <- c("Quadrant", "N (Resection/LRST)", "IPTW 5-yr RD (95% CI)", "Outcome-Adjusted RD", "Adjusted HR (95% CI)", "P Value")

cat("\n=== Supplementary Table S6 ===\n")
print(s6_table, row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. IPTW-Adjusted Kaplan–Meier Survival Curves (Figure 2A & 2B)
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Figure 2A (OS) and Figure 2B (RFS)]

cat("\n══ Section 9: Generating Figure 2A (OS) & Figure 2B (RFS) ══\n")

# Auto-detect analysis dataset (dat_pooled or dat_scored)
dat_km_base <- if (exists("dat_pooled")) dat_pooled else dat_scored

# Figure 2A: IPTW-Adjusted Overall Survival (OS) Curves across Quadrants

# 1. OS Data Preparation
dat_km_prep_os <- dat_km_base %>%
  filter(!is.na(quadrant) & !is.na(w_ATT_win) & !is.na(OS_time_5y)) %>%
  mutate(time_yr = if (max(OS_time_5y, na.rm = TRUE) > 10) OS_time_5y / 365.25 else OS_time_5y)

# 2. OS KM Plot Generator Function
make_single_km_os_plot <- function(df, quad_key, panel_title, show_xlab = TRUE, show_ylab = TRUE) {
  sub_df <- df %>% filter(quadrant == quad_key)
  n_r <- sum(sub_df$trt == 1, na.rm = TRUE)
  n_l <- sum(sub_df$trt == 0, na.rm = TRUE)
  
  km_r <- survfit(Surv(time_yr, OS_status_5y) ~ 1, weights = w_ATT_win, data = filter(sub_df, trt == 1))
  km_l <- survfit(Surv(time_yr, OS_status_5y) ~ 1, weights = w_ATT_win, data = filter(sub_df, trt == 0))
  
  lbl_r <- sprintf("Resection (n=%d)", n_r)
  lbl_l <- sprintf("LRST (n=%d)", n_l)
  
  df_plot <- bind_rows(
    data.frame(time = c(0, km_r$time), surv = c(1, km_r$surv) * 100, trt_group = lbl_r),
    data.frame(time = c(0, km_l$time), surv = c(1, km_l$surv) * 100, trt_group = lbl_l)
  ) %>% mutate(trt_group = factor(trt_group, levels = c(lbl_r, lbl_l)))
  
  stat_annot <- tryCatch({
    s_cox <- summary(coxph(Surv(time_yr, OS_status_5y) ~ trt, data = sub_df, weights = w_ATT_win, robust = TRUE))
    p_val <- s_cox$coefficients[1, "Pr(>|z|)"]
    p_str <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
    sprintf("IPTW HR: %.2f (%.2f–%.2f)\n%s", s_cox$coefficients[1, "exp(coef)"], s_cox$conf.int[1, "lower .95"], s_cox$conf.int[1, "upper .95"], p_str)
  }, error = function(e) "")
  
  ggplot(df_plot, aes(x = time, y = surv, color = trt_group)) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "#a1a1a1", linewidth = 0.5) +
    geom_step(linewidth = 0.95, alpha = 0.95) +
    annotate("text", x = 4.85, y = 18, label = stat_annot, hjust = 1, vjust = 0, size = 3.2, fontface = "bold", color = "#222222") +
    scale_color_manual(values = setNames(c("#1a5fb4", "#c01c28"), c(lbl_r, lbl_l))) +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5, expand = c(0.005, 0.005)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20), expand = c(0.005, 0.005)) +
    labs(
      title = panel_title,
      x     = if (show_xlab) "Years from first HCC-directed therapy" else NULL,
      y     = if (show_ylab) "IPTW-adjusted overall survival (%)" else NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(size = 10.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 6)),
      panel.grid.major = element_line(color = "#f2f2f2", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.text        = element_text(size = 9.5, color = "black"),
      axis.title       = element_text(size = 10.5, color = "black"),
      axis.title.y     = element_text(margin = margin(r = 6)),
      axis.title.x     = element_text(margin = margin(t = 6)),
      legend.position  = c(0.03, 0.05),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.width  = unit(1.2, "line"),
      legend.text       = element_text(size = 9.2, face = "bold", color = "#222222"),
      legend.spacing.y  = unit(0.08, "cm"),
      legend.title      = element_blank()
    )
}

# 3. Assemble Figure 2A (OS) Panels
p_os_tl <- make_single_km_os_plot(dat_km_prep_os, "Low-O / High-S", "Low O / High S (surgical futility)", FALSE, TRUE)
p_os_tr <- make_single_km_os_plot(dat_km_prep_os, "High-O / High-S", "High O / High S (doomed)", FALSE, FALSE)
p_os_bl <- make_single_km_os_plot(dat_km_prep_os, "Low-O / Low-S", "Low O / Low S (standard resection)", TRUE, TRUE)
p_os_br <- make_single_km_os_plot(dat_km_prep_os, "High-O / Low-S", "High O / Low S (oncological futility)", TRUE, FALSE)

p_fig2a_os <- (p_os_tl | p_os_tr) / (p_os_bl | p_os_br) +
  plot_annotation(
    title = "Figure 2A. IPTW-adjusted overall survival by O/S quadrant (resection vs LRST)",
    theme = theme(plot.title = element_text(size = 13.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 10, t = 5)))
  )

print(p_fig2a_os)


# Figure 2B: IPTW-Adjusted Recurrence-Free Survival (RFS) Curves across Quadrants

# 1. RFS Data Preparation
dat_km_prep_rfs <- dat_km_base %>%
  filter(!is.na(quadrant) & !is.na(w_ATT_win) & !is.na(RFS_time_5y)) %>%
  mutate(time_yr = if (max(RFS_time_5y, na.rm = TRUE) > 10) RFS_time_5y / 365.25 else RFS_time_5y)

# 2. RFS KM Plot Generator Function
make_single_km_rfs_plot <- function(df, quad_key, panel_title, show_xlab = TRUE, show_ylab = TRUE) {
  sub_df <- df %>% filter(quadrant == quad_key)
  n_r <- sum(sub_df$trt == 1, na.rm = TRUE)
  n_l <- sum(sub_df$trt == 0, na.rm = TRUE)
  
  km_r <- survfit(Surv(time_yr, RFS_status_5y) ~ 1, weights = w_ATT_win, data = filter(sub_df, trt == 1))
  km_l <- survfit(Surv(time_yr, RFS_status_5y) ~ 1, weights = w_ATT_win, data = filter(sub_df, trt == 0))
  
  lbl_r <- sprintf("Resection (n=%d)", n_r)
  lbl_l <- sprintf("LRST (n=%d)", n_l)
  
  df_plot <- bind_rows(
    data.frame(time = c(0, km_r$time), surv = c(1, km_r$surv) * 100, trt_group = lbl_r),
    data.frame(time = c(0, km_l$time), surv = c(1, km_l$surv) * 100, trt_group = lbl_l)
  ) %>% mutate(trt_group = factor(trt_group, levels = c(lbl_r, lbl_l)))
  
  stat_annot <- tryCatch({
    s_cox <- summary(coxph(Surv(time_yr, RFS_status_5y) ~ trt, data = sub_df, weights = w_ATT_win, robust = TRUE))
    p_val <- s_cox$coefficients[1, "Pr(>|z|)"]
    p_str <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
    sprintf("IPTW HR: %.2f (%.2f–%.2f)\n%s", s_cox$coefficients[1, "exp(coef)"], s_cox$conf.int[1, "lower .95"], s_cox$conf.int[1, "upper .95"], p_str)
  }, error = function(e) "")
  
  ggplot(df_plot, aes(x = time, y = surv, color = trt_group)) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "#a1a1a1", linewidth = 0.5) +
    geom_step(linewidth = 0.95, alpha = 0.95) +
    annotate("text", x = 4.85, y = 18, label = stat_annot, hjust = 1, vjust = 0, size = 3.2, fontface = "bold", color = "#222222") +
    scale_color_manual(values = setNames(c("#1a5fb4", "#c01c28"), c(lbl_r, lbl_l))) +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5, expand = c(0.005, 0.005)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20), expand = c(0.005, 0.005)) +
    labs(
      title = panel_title,
      x     = if (show_xlab) "Years from first HCC-directed therapy" else NULL,
      y     = if (show_ylab) "IPTW-adjusted recurrence-free survival (%)" else NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(size = 10.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 6)),
      panel.grid.major = element_line(color = "#f2f2f2", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.text        = element_text(size = 9.5, color = "black"),
      axis.title       = element_text(size = 10.5, color = "black"),
      axis.title.y     = element_text(margin = margin(r = 6)),
      axis.title.x     = element_text(margin = margin(t = 6)),
      legend.position  = c(0.03, 0.05),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.width  = unit(1.2, "line"),
      legend.text       = element_text(size = 9.2, face = "bold", color = "#222222"),
      legend.spacing.y  = unit(0.08, "cm"),
      legend.title      = element_blank()
    )
}

# 3. Assemble Figure 2B (RFS) Panels
p_rfs_tl <- make_single_km_rfs_plot(dat_km_prep_rfs, "Low-O / High-S", "Low O / High S (surgical futility)", FALSE, TRUE)
p_rfs_tr <- make_single_km_rfs_plot(dat_km_prep_rfs, "High-O / High-S", "High O / High S (doomed)", FALSE, FALSE)
p_rfs_bl <- make_single_km_rfs_plot(dat_km_prep_rfs, "Low-O / Low-S", "Low O / Low S (standard resection)", TRUE, TRUE)
p_rfs_br <- make_single_km_rfs_plot(dat_km_prep_rfs, "High-O / Low-S", "High O / Low S (oncological futility)", TRUE, FALSE)

p_fig2b_rfs <- (p_rfs_tl | p_rfs_tr) / (p_rfs_bl | p_rfs_br) +
  plot_annotation(
    title = "Figure 2B. IPTW-adjusted recurrence-free survival by O/S quadrant (resection vs LRST)",
    theme = theme(plot.title = element_text(size = 13.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 10, t = 5)))
  )

print(p_fig2b_rfs)

# ------------------------------------------------------------------------------
# 10. Table 3: Primary Counterfactual Outcomes Assembly
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Table 3]
phenotype_map <- c(
  "Low-O / Low-S"   = "Standard resection",
  "High-O / Low-S"  = "Oncological futility",
  "Low-O / High-S"  = "Surgical futility",
  "High-O / High-S" = "Dual elevation"
)

table3_observed <- quad_summary_final %>%
  mutate(
    Quadrant                      = as.character(`Quadrant Name`),
    Phenotype                     = unname(phenotype_map[as.character(`Quadrant Name`)]),
    `n (% of cohort)`             = sprintf("%s (%s)", format(Count_N, big.mark = ","), Pct),
    `Observed 5-yr OS`            = gsub(" .*", "", `5-yr OS (95% CI)`),
    `5-yr CIF recurrence`         = `5-yr Recurrence CIF`,
    `5-yr CIF non-recurrent death` = `5-yr Non-recurrent Death CIF`
  ) %>%
  select(
    Quadrant, Phenotype, `n (% of cohort)`, 
    `Observed 5-yr OS`, `5-yr CIF recurrence`, `5-yr CIF non-recurrent death`
  )

table3_att <- fig1_results %>%
  mutate(
    Quadrant                      = as.character(quadrant),
    `LRST Comparator, n (ESS)`    = sprintf("%d (%.1f)", n_lrst, ess_lrst),
    `ATT 5-yr RD vs LRST (95% CI)` = sprintf("%+.1f (%+.1f, %+.1f)", OS_diff_pp, CI_lo, CI_hi)
  ) %>%
  select(Quadrant, `LRST Comparator, n (ESS)`, `ATT 5-yr RD vs LRST (95% CI)`)
      
table3_final <- table3_observed %>% left_join(table3_att, by = "Quadrant")

cat("\n=== Table 3. Main Cohort Outcomes and ATT 5-Year Risk Differences ===\n")
print(table3_final, row.names = FALSE)

# ------------------------------------------------------------------------------
# 11. Supplementary Table S3: Multi-Year Bootstrap Treatment Effect Estimates
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Table S3]
generate_supp_table_s3_standard <- function(df, ps_formula, n_boot = 1000, seed = 2026) {
  
  set.seed(seed)
  horizons_days  <- c(2 * 365.25, 3 * 365.25, 5 * 365.25) # 2y, 3y, 5y in days
  horizon_labels <- c(2, 3, 5)
  quads          <- levels(df$quadrant)
  
  # Helper function: Calculate point RD at time t
  calc_quad_rd_at_t <- function(sub_df, t_days) {
    get_os <- function(trt_val) {
      d <- sub_df %>% filter(trt == trt_val)
      if (nrow(d) < 2) return(NA_real_)
      
      km <- tryCatch(
        survfit(Surv(OS_time_5y, OS_status_5y) ~ 1, weights = d$w_ATT_win, data = d),
        error = function(e) NULL
      )
      if (is.null(km)) return(NA_real_)
      
      s_sum <- summary(km, times = t_days, extend = TRUE)
      if (length(s_sum$surv) > 0) s_sum$surv[1] * 100 else NA_real_
    }
    
    os_r <- get_os(1)
    os_l <- get_os(0)
    if (!is.na(os_r) && !is.na(os_l)) os_r - os_l else NA_real_
  }
  
  # 1. Point Estimates on Original Data
  point_list <- list()
  for (q in quads) {
    sub_q <- df %>% filter(quadrant == q & !is.na(w_ATT_win))
    for (i in seq_along(horizons_days)) {
      rd_val <- calc_quad_rd_at_t(sub_q, horizons_days[i])
      point_list[[paste(q, horizon_labels[i], sep = "_")]] <- rd_val
    }
  }
  
  # 2. Stratified Bootstrap with PS Model Re-fitting (1000 Replicates)
  cat(sprintf("── Running Stratified Bootstrap-%d with PS model re-fitting...\n", n_boot))
  
  matrix_cols <- names(point_list)
  boot_mat    <- matrix(NA, nrow = n_boot, ncol = length(matrix_cols))
  colnames(boot_mat) <- matrix_cols
  
  for (b in seq_len(n_boot)) {
    # Stratified resample by quadrant x treatment
    boot_df <- df %>%
      group_by(quadrant, trt) %>%
      slice_sample(prop = 1, replace = TRUE) %>%
      ungroup()
    
    # Re-fit Propensity Score model and update ATT weights
    ps_b <- tryCatch({
      suppressWarnings({
        m_b <- glm(ps_formula, data = boot_df, family = binomial(link = "logit"))
        predict(m_b, newdata = boot_df, type = "response")
      })
    }, error = function(e) boot_df$ps)
    
    ps_b_clean        <- pmin(pmax(ps_b, 0.01), 0.99)
    w_ATT_b_raw       <- ifelse(boot_df$trt == 1, 1.0, ps_b_clean / (1 - ps_b_clean))
    boot_df$w_ATT_win <- pmin(w_ATT_b_raw, quantile(w_ATT_b_raw, 0.99, na.rm = TRUE))
    
    # Calculate RD for each quadrant and horizon in bootstrap sample
    for (q in quads) {
      sub_q_b <- boot_df %>% filter(quadrant == q)
      for (i in seq_along(horizons_days)) {
        col_key <- paste(q, horizon_labels[i], sep = "_")
        boot_mat[b, col_key] <- calc_quad_rd_at_t(sub_q_b, horizons_days[i])
      }
    }
    
    if (b %% 200 == 0) cat(sprintf("   Bootstrap %d/%d completed\n", b, n_boot))
  }
  
  # 3. Assemble Supplementary Table S3
  results_list <- list()
  for (q in quads) {
    for (i in seq_along(horizons_days)) {
      col_key <- paste(q, horizon_labels[i], sep = "_")
      
      rd_point  <- point_list[[col_key]]
      boot_vals <- boot_mat[, col_key]
      boot_vals <- boot_vals[!is.na(boot_vals)]
      
      ci_lo <- quantile(boot_vals, probs = 0.025, na.rm = TRUE)
      ci_hi <- quantile(boot_vals, probs = 0.975, na.rm = TRUE)
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Quadrant                      = gsub(" / ", "/", q),
        `Horizon (Years)`            = horizon_labels[i],
        `ATT RD (Percentage Points)` = sprintf("%+.1f", rd_point),
        `95% Bootstrap CI`           = sprintf("%+.1f to %+.1f", ci_lo, ci_hi),
        `n_boot`                     = n_boot,
        check.names                  = FALSE,
        stringsAsFactors             = FALSE
      )
    }
  }
  
  bind_rows(results_list)
}

# Run and format output
dat_supp_base <- if (exists("dat_pooled")) dat_pooled else dat_scored

supp_table_s3_final <- generate_supp_table_s3_standard(
  df         = dat_supp_base, 
  ps_formula = ps_formula, 
  n_boot     = 1000, 
  seed       = 2026
)

cat("\n=== Supplementary Table S3. Multi-Year Bootstrap Treatment Effect Estimates ===\n\n")
print(supp_table_s3_final, row.names = FALSE)

# ------------------------------------------------------------------------------
# 12. Supplementary Tables S4 & S5: Sensitivity Analyses
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Tables S4 and S5]
generate_supp_table_s4_trim <- function(df) {
  t_2y <- 730.5; t_3y <- 1095.75; t_5y <- 1826.25
  get_os_pct <- function(sub_data, time_pt) {
    if (nrow(sub_data) < 2) return(NA_real_)
    km <- tryCatch(survfit(Surv(OS_time_5y, OS_status_5y) ~ 1, weights = w_ATT_win, data = sub_data), error = function(e) NULL)
    if (is.null(km)) return(NA_real_)
    s_sum <- summary(km, times = time_pt, extend = TRUE)
    if (length(s_sum$surv) > 0) return(s_sum$surv[1] * 100) else return(NA_real_)
  }
  df_trim <- df %>% filter(ps >= 0.10 & ps <= 0.90)
  table_rows <- list()
  
  for (q in levels(df$quadrant)) {
    sub_trim <- df_trim %>% filter(quadrant == q & !is.na(w_ATT_win) & !is.na(OS_time_5y))
    d_res_t  <- sub_trim %>% filter(trt == 1); d_lrst_t <- sub_trim %>% filter(trt == 0)
    
    os_2y_r <- get_os_pct(d_res_t, t_2y); os_2y_l <- get_os_pct(d_lrst_t, t_2y)
    os_3y_r <- get_os_pct(d_res_t, t_3y); os_3y_l <- get_os_pct(d_lrst_t, t_3y)
    os_5y_r <- get_os_pct(d_res_t, t_5y); os_5y_l <- get_os_pct(d_lrst_t, t_5y)
    
    fmt_os_pair <- function(r, l) { if (is.na(r) || is.na(l)) "N/A" else sprintf("%.1f / %.1f", r, l) }
    fmt_rd      <- function(r, l) { if (is.na(r) || is.na(l)) "N/A" else sprintf("%.1f", r - l) }
    
    table_rows[[q]] <- data.frame(
      Quadrant = gsub(" / ", "/", q), `Resection (n)` = nrow(d_res_t), `LRST (n)` = nrow(d_lrst_t), `LRST ESS` = calc_ess(d_lrst_t$w_ATT_win),
      `2-y OS (%) R/L` = fmt_os_pair(os_2y_r, os_2y_l), `RD (pp)_2y` = fmt_rd(os_2y_r, os_2y_l),
      `3-y OS (%) R/L` = fmt_os_pair(os_3y_r, os_3y_l), `RD (pp)_3y` = fmt_rd(os_3y_r, os_3y_l),
      `5-y OS (%) R/L` = fmt_os_pair(os_5y_r, os_5y_l), `RD (pp)_5y` = fmt_rd(os_5y_r, os_5y_l),
      check.names = FALSE
    )
  }
  bind_rows(table_rows)
}

supp_table_s4_trim_expert <- generate_supp_table_s4_trim(dat_pooled)
cat("\n=== Supplementary Table S4 ===\n")
print(supp_table_s4_trim_expert, row.names = FALSE)

generate_supp_table_s5_tace <- function(df, ps_covars) {
  dat_tace_subset <- df %>% filter(trt == 1 | (trt == 0 & (if ("TACE" %in% colnames(.)) TACE == 1 else LRST_type == "TACE")))
  ps_formula <- as.formula(paste("trt ~", paste(ps_covars, collapse = " + ")))
  ps_model_tace <- glm(ps_formula, data = dat_tace_subset, family = binomial(link = "logit"))
  dat_tace_subset$ps_tace <- predict(ps_model_tace, newdata = dat_tace_subset, type = "response")
  
  dat_tace_subset <- dat_tace_subset %>% mutate(w_att_tace_raw = ifelse(trt == 1, 1, ps_tace / (1 - ps_tace)))
  w_cap <- quantile(dat_tace_subset$w_att_tace_raw[dat_tace_subset$trt == 0], probs = 0.99, na.rm = TRUE)
  dat_tace_subset <- dat_tace_subset %>% mutate(w_att_tace_win = ifelse(trt == 1, 1, pmin(w_att_tace_raw, w_cap)))
  
  t_2y <- 730.5; t_3y <- 1095.75; t_5y <- 1826.25
  table_rows <- list()
  
  for (q in levels(df$quadrant)) {
    sub_q  <- dat_tace_subset %>% filter(quadrant == q & !is.na(w_att_tace_win) & !is.na(OS_time_5y))
    d_res  <- sub_q %>% filter(trt == 1); d_tace <- sub_q %>% filter(trt == 0)
    
    get_os_pct <- function(data_sub, time_pt) {
      if (nrow(data_sub) < 2) return(NA_real_)
      km <- tryCatch(survfit(Surv(OS_time_5y, OS_status_5y) ~ 1, weights = w_att_tace_win, data = data_sub), error = function(e) NULL)
      if (is.null(km)) return(NA_real_)
      s_sum <- summary(km, times = time_pt, extend = TRUE)
      if (length(s_sum$surv) > 0) return(s_sum$surv[1] * 100) else return(NA_real_)
    }
    
    os_2y_r <- get_os_pct(d_res, t_2y); os_2y_t <- get_os_pct(d_tace, t_2y)
    os_3y_r <- get_os_pct(d_res, t_3y); os_3y_t <- get_os_pct(d_tace, t_3y)
    os_5y_r <- get_os_pct(d_res, t_5y); os_5y_t <- get_os_pct(d_tace, t_5y)
    
    fmt_os_pair <- function(r, t) { if (is.na(r) || is.na(t)) "N/A" else sprintf("%.1f / %.1f", r, t) }
    fmt_rd      <- function(r, t) { if (is.na(r) || is.na(t)) "N/A" else sprintf("%+.1f", r - t) }
    
    table_rows[[q]] <- data.frame(
      `Risk phenotype` = gsub(" / ", "/", q), `Resection n` = nrow(d_res), `TACE n` = nrow(d_tace), `TACE ESS` = calc_ess(d_tace$w_att_tace_win),
      `2-y OS (%) R/T` = fmt_os_pair(os_2y_r, os_2y_t), `2-y RD (pp)` = fmt_rd(os_2y_r, os_2y_t),
      `3-y OS (%) R/T` = fmt_os_pair(os_3y_r, os_3y_t), `3-y RD (pp)` = fmt_rd(os_3y_r, os_3y_t),
      `5-y OS (%) R/T` = fmt_os_pair(os_5y_r, os_5y_t), `5-y RD (pp)` = fmt_rd(os_5y_r, os_5y_t),
      check.names = FALSE
    )
  }
  bind_rows(table_rows)
}

supp_table_s5_tace_full <- generate_supp_table_s5_tace(dat_pooled, ps_vars_in_data)
cat("\n=== Supplementary Table S5 ===\n")
print(supp_table_s5_tace_full, row.names = FALSE)

# ------------------------------------------------------------------------------
# 13. Supplementary Figure S4: Forest Plot (Doubly Adjusted Counterfactual Effects)
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Figure S4]
df_single <- data.frame(
  id            = 4:1,
  Quadrant      = c("Low-O / Low-S", "High-O / Low-S", "Low-O / High-S", "High-O / High-S"),
  Clinical_Note = c("Standard resection", "Oncological futility", "Surgical futility", "Dual futility"),
  N_Label       = c("1804 / 218", "786 / 166", "329 / 211", "128 / 116"),
  RD_Label      = c("+20.5 pp", "+14.8 pp", "+18.4 pp", "+15.8 pp"),
  HR            = c(0.30, 0.45, 0.52, 0.43),
  HR_lo         = c(0.18, 0.25, 0.31, 0.17),
  HR_hi         = c(0.51, 0.82, 0.87, 1.08),
  HR_CI_Text    = c("0.30 (0.18–0.51)", "0.45 (0.25–0.82)", "0.52 (0.31–0.87)", "0.43 (0.17–1.08)"),
  P_val_Text    = c("< 0.001", "0.010", "0.012", "0.074"),
  Is_Sig        = c(TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

p_top_journal <- ggplot(df_single, aes(y = id)) +
  geom_rect(aes(ymin = id - 0.5, ymax = id + 0.5, xmin = -2.12, xmax = 2.12), fill = rep(c("white", "#F8F9FA"), length.out = 4), alpha = 0.9) +
  geom_text(aes(x = -2.10, label = Quadrant), fontface = "bold", hjust = 0, size = 3.6) +
  geom_text(aes(x = -1.45, label = Clinical_Note), fontface = "italic", color = "grey30", hjust = 0, size = 3.2) +
  geom_text(aes(x = -0.85, label = N_Label), hjust = 0.5, size = 3.3) +
  geom_text(aes(x = -0.30, label = RD_Label), hjust = 0.5, size = 3.3) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey35", size = 0.55) +
  geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi, color = Is_Sig), height = 0.18, size = 0.8) +
  geom_point(aes(x = HR, color = Is_Sig, fill = Is_Sig), shape = 22, size = 3.5, stroke = 0.6) +
  geom_text(aes(x = 1.55, label = HR_CI_Text, fontface = ifelse(Is_Sig, "bold", "plain")), hjust = 0.5, size = 3.3) +
  geom_text(aes(x = 2.05, label = P_val_Text, fontface = ifelse(Is_Sig, "bold", "plain")), hjust = 0.5, size = 3.3) +
  scale_color_manual(values = c("TRUE" = "#00A087", "FALSE" = "#E64B35")) +
  scale_fill_manual(values = c("TRUE" = "#00A087", "FALSE" = "#E64B35")) +
  annotate("text", x = c(-2.10, -1.45, -0.85, -0.30, 0.65, 1.55, 2.05), y = 4.8,
           label = c("Risk Quadrant", "Clinical Profile", "N (Resect/LRST)", "Outcome-Adj RD", "Adjusted HR (95% CI)", "Adjusted HR (95% CI)", "P Value"),
           fontface = "bold", size = 3.6, hjust = c(0, 0, 0.5, 0.5, 0.5, 0.5, 0.5)) +
  annotate("segment", x = -2.10, xend = 2.12, y = 5.2, yend = 5.2, size = 0.8, color = "black") +
  annotate("segment", x = -2.10, xend = 2.12, y = 4.35, yend = 4.35, size = 0.6, color = "black") +
  annotate("segment", x = -2.10, xend = 2.12, y = 0.55, yend = 0.55, size = 0.8, color = "black") +
  annotate("segment", x = 0.2, xend = 1.2, y = 0.3, yend = 0.3, color = "black", size = 0.5) +
  annotate("segment", x = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), xend = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), y = 0.3, yend = 0.22, color = "black", size = 0.5) +
  annotate("text", x = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), y = 0.05, label = c("0.2", "0.4", "0.6", "0.8", "1.0", "1.2"), size = 2.8, fontface = "bold") +
  annotate("text", x = 0.65, y = -0.22, label = "← Favors Resection | Favors LRST →", size = 3.0, fontface = "bold", color = "grey20") +
  coord_cartesian(xlim = c(-2.15, 2.15), ylim = c(-0.35, 5.5), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 15, r = 20, b = 15, l = 20), legend.position = "none")

final_fig_pretty <- p_top_journal +
  labs(
    title = "Supplementary Figure 4. Counterfactual Treatment Effects Across Risk Quadrants (Doubly Adjusted)",
    subtitle = "IPTW-weighted multivariable outcome Cox regression adjusting for residual imbalanced baseline covariates",
    caption = "Treatment-by-quadrant interaction: P = 0.027. Outcome models explicitly adjusted for residual imbalanced covariates (|SMD| > 0.10: tumor diameter, HBsAg status, cirrhosis, and serum AFP).\nGreen squares denote statistically significant survival benefit (P < 0.05, 95% CI excludes 1.0); red square denotes surgical futility (P = 0.074, 95% CI includes 1.0)."
  ) +
  theme(plot.title = element_text(face = "bold", size = 12.5, color = "black"),
        plot.subtitle = element_text(size = 9.5, color = "grey30"),
        plot.caption = element_text(size = 8.5, color = "black", face = "italic", hjust = 0))

print(final_fig_pretty)
