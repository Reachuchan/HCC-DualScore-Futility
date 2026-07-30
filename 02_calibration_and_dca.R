# ==============================================================================
# Script Name: 02_calibration_and_dca.R
# Purpose: Model calibration curves and Decision Curve Analysis (DCA).
# Outputs: Supplementary Figure S1, Supplementary Figure S2.
# Language: R (>= 4.0.0)
# Requirements: Requires workspace object 'dat_scored' from Script 01.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(survival)
library(cmprsk)
library(patchwork)
library(dcurves)
library(tidyverse)

# ------------------------------------------------------------------------------
# 2. Supplementary Figure S1: Model Calibration Curves
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Figure S1]
get_calibration_data <- function(data, x_vars, time_var, event_var, cause_code) {
  data_clean <- data %>% filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]), .data[[time_var]] > 0)
  
  lp <- get_lp_ridge(data_clean, data_clean, x_vars, time_var, event_var, lambda_val = 0.5)$lp
  fit_cox <- coxph(Surv(data_clean[[time_var]], data_clean[[event_var]]) ~ lp, data = data_clean)
  
  t_5yr <- ifelse(max(data_clean[[time_var]], na.rm = TRUE) > 100, 1826.25, 5)
  pred_cif <- (1 - as.numeric(summary(survfit(fit_cox, newdata = data_clean), times = t_5yr)$surv)) * 100
  
  data_cal <- data.frame(pred = pred_cif, time = as.numeric(data_clean[[time_var]]), status = as.numeric(data_clean[[event_var]])) %>%
    filter(!is.na(pred), !is.na(time), !is.na(status)) %>%
    mutate(decile = ntile(pred, 10))
  
  res <- data_cal %>%
    group_by(decile) %>%
    do({
      d_sub <- .
      mean_pred <- mean(d_sub$pred, na.rm = TRUE)
      
      obs_cif <- tryCatch({
        ci_fit <- cuminc(ftime = d_sub$time, fstatus = d_sub$status, cencode = 0)
        key_name <- paste("1", cause_code)
        if (key_name %in% names(ci_fit)) {
          t_vals <- ci_fit[[key_name]]$time
          idx <- max(which(t_vals <= t_5yr))
          if (is.finite(idx) && idx > 0) ci_fit[[key_name]]$est[idx] * 100 else 0
        } else if (length(ci_fit) > 0) {
          first_key <- names(ci_fit)[1]
          t_vals <- ci_fit[[first_key]]$time
          idx <- max(which(t_vals <= t_5yr))
          if (is.finite(idx) && idx > 0) ci_fit[[first_key]]$est[idx] * 100 else 0
        } else { 0 }
      }, error = function(e) {
        mean(d_sub$status == cause_code & d_sub$time <= t_5yr, na.rm = TRUE) * 100
      })
      data.frame(mean_pred = mean_pred, obs_cif = obs_cif)
    }) %>%
    ungroup()
  return(res)
}

cal_o <- get_calibration_data(dat_scored, O_VARS, "O_time", "O_event", cause_code = 1)
cal_s <- get_calibration_data(dat_scored, S_VARS, "S_time", "S_event", cause_code = 2)

theme_target <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major  = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title        = element_text(size = 10, face = "plain", hjust = 0.5, margin = margin(b = 8)),
    axis.title        = element_text(size = 10, face = "plain"),
    axis.text         = element_text(size = 9, color = "black"),
    legend.title      = element_blank(),
    legend.position   = c(0.24, 0.92),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

ref_line_df <- data.frame(x = c(0, 80), y = c(0, 80), type = "Perfect calibration")

p_s1a <- ggplot() +
  geom_line(data = ref_line_df, aes(x = x, y = y, linetype = type), color = "gray40", linewidth = 0.7) +
  geom_point(data = cal_o, aes(x = mean_pred, y = obs_cif), shape = 21, color = "#c0262b", fill = "white", stroke = 1.3, size = 2.5) +
  scale_linetype_manual(values = c("Perfect calibration" = "dashed")) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 10), expand = c(0, 0)) +
  labs(title = "O-score: predicted vs observed 5-y CIF for first recurrence", x = "Mean predicted 5-year CIF (%)", y = "Observed 5-year CIF (%)") +
  theme_target

p_s1b <- ggplot() +
  geom_line(data = ref_line_df, aes(x = x, y = y, linetype = type), color = "gray40", linewidth = 0.7) +
  geom_point(data = cal_s, aes(x = mean_pred, y = obs_cif), shape = 21, color = "#0e58c7", fill = "white", stroke = 1.3, size = 2.5) +
  scale_linetype_manual(values = c("Perfect calibration" = "dashed")) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 10), expand = c(0, 0)) +
  labs(title = "S-score: predicted vs observed 5-y CIF for non-recurrent death", x = "Mean predicted 5-year CIF (%)", y = "Observed 5-year CIF (%)") +
  theme_target

fig_s1 <- (p_s1a | p_s1b) + plot_annotation(title = "Calibration of dual-score models — decile of predicted 5-year CIF vs observed", theme = theme(plot.title = element_text(size = 12, face = "plain", hjust = 0.5)))
print(fig_s1)
# ggsave("Supplementary_Figure_S1.pdf", fig_s1, width = 10, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 3. Supplementary Figure S2: Vickers Decision Curve Analysis
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Figure S2]
calc_vickers_dca <- function(data, x_vars, time_var, event_var, cause_code, thresholds = seq(0.05, 0.85, by = 0.02)) {
  data_clean <- data %>% filter(!is.na(.data[[time_var]]), .data[[time_var]] > 0, !is.na(.data[[event_var]]))
  
  lp <- get_lp_ridge(data_clean, data_clean, x_vars, time_var, event_var, lambda_val = 0.5)$lp
  fit_cox <- coxph(Surv(data_clean[[time_var]], data_clean[[event_var]]) ~ lp, data = data_clean)
  
  t_5yr <- ifelse(max(data_clean[[time_var]], na.rm = TRUE) > 100, 1826.25, 5)
  pred_prob <- 1 - as.numeric(summary(survfit(fit_cox, newdata = data_clean), times = t_5yr)$surv)
  N <- nrow(data_clean)
  
  p_overall <- tryCatch({
    ci_all <- cuminc(ftime = data_clean[[time_var]], fstatus = data_clean[[event_var]], cencode = 0)
    key_name <- paste("1", cause_code)
    k_use <- if (key_name %in% names(ci_all)) key_name else names(ci_all)[1]
    idx <- max(which(ci_all[[k_use]]$time <= t_5yr))
    ci_all[[k_use]]$est[idx]
  }, error = function(e) mean(data_clean[[event_var]] == cause_code & data_clean[[time_var]] <= t_5yr))
  
  dca_res <- data.frame()
  for (pt in thresholds) {
    flagged <- pred_prob >= pt
    n_flagged <- sum(flagged, na.rm = TRUE)
    if (n_flagged > 0) {
      sub_data <- data_clean[flagged, ]
      p_sub <- tryCatch({
        ci_sub <- cuminc(ftime = sub_data[[time_var]], fstatus = sub_data[[event_var]], cencode = 0)
        key_name <- paste("1", cause_code)
        k_use <- if (key_name %in% names(ci_sub)) key_name else names(ci_sub)[1]
        idx <- max(which(ci_sub[[k_use]]$time <= t_5yr))
        ci_sub[[k_use]]$est[idx]
      }, error = function(e) mean(sub_data[[event_var]] == cause_code & sub_data[[time_var]] <= t_5yr))
      tp <- p_sub * (n_flagged / N)
      fp <- (1 - p_sub) * (n_flagged / N)
      nb_model <- tp - fp * (pt / (1 - pt))
    } else { nb_model <- 0 }
    nb_all <- p_overall - (1 - p_overall) * (pt / (1 - pt))
    nb_none <- 0
    dca_res <- rbind(dca_res, data.frame(pt = pt * 100, nb_model = nb_model, nb_all = nb_all, nb_none = nb_none))
  }
  return(dca_res)
}

dca_o <- calc_vickers_dca(dat_scored, O_VARS, "O_time", "O_event", cause_code = 1)
dca_s <- calc_vickers_dca(dat_scored, S_VARS, "S_time", "S_event", cause_code = 2)

legend_levels <- c("Score-based decision", "Treat all (flag every patient)", "Treat none (flag no patient)")

dca_o_long <- bind_rows(
  dca_o %>% select(pt, nb = nb_model) %>% mutate(type = "Score-based decision"),
  dca_o %>% select(pt, nb = nb_all)   %>% mutate(type = "Treat all (flag every patient)"),
  dca_o %>% select(pt, nb = nb_none)  %>% mutate(type = "Treat none (flag no patient)")
) %>% mutate(type = factor(type, levels = legend_levels))

dca_s_long <- bind_rows(
  dca_s %>% select(pt, nb = nb_model) %>% mutate(type = "Score-based decision"),
  dca_s %>% select(pt, nb = nb_all)   %>% mutate(type = "Treat all (flag every patient)"),
  dca_s %>% select(pt, nb = nb_none)  %>% mutate(type = "Treat none (flag no patient)")
) %>% mutate(type = factor(type, levels = legend_levels))

theme_dca_target <- theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(color = "gray92", linewidth = 0.3), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), legend.position = c(0.70, 0.83))

p_s2a <- ggplot(dca_o_long, aes(x = pt, y = nb, color = type, linetype = type)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("Score-based decision" = "#c0262b", "Treat all (flag every patient)" = "gray50", "Treat none (flag no patient)" = "black")) +
  scale_linetype_manual(values = c("Score-based decision" = "solid", "Treat all (flag every patient)" = "dashed", "Treat none (flag no patient)" = "dotted")) +
  scale_x_continuous(limits = c(5, 90), breaks = seq(10, 90, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.05, 0.5), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
  labs(title = "O-score: net benefit at varying 5-y recurrence-risk thresholds", x = "Threshold probability of 5-year cause-specific event (%)", y = "Net benefit") +
  theme_dca_target

p_s2b <- ggplot(dca_s_long, aes(x = pt, y = nb, color = type, linetype = type)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("Score-based decision" = "#0e58c7", "Treat all (flag every patient)" = "gray50", "Treat none (flag no patient)" = "black")) +
  scale_linetype_manual(values = c("Score-based decision" = "solid", "Treat all (flag every patient)" = "dashed", "Treat none (flag no patient)" = "dotted")) +
  scale_x_continuous(limits = c(5, 90), breaks = seq(10, 90, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.05, 0.5), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
  labs(title = "S-score: net benefit at varying 5-y non-recurrent-mortality thresholds", x = "Threshold probability of 5-year cause-specific event (%)", y = "Net benefit") +
  theme_dca_target

fig_s2 <- (p_s2a | p_s2b) + plot_annotation(title = "Decision-curve analysis (Vickers) — dual-score framework")
print(fig_s2)
# ggsave("Supplementary_Figure_S2.pdf", fig_s2, width = 10, height = 5, dpi = 300)
