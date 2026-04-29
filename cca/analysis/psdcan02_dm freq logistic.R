rm(list=ls());gc();source(".Rprofile")

analytic_df_wide <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/archive/psdcpre04_analytic dataset wide.RDS"))

analytic_df_wide <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre05_analytic dataset wide.RDS"))

# ============================================================================
# Frequencies and simple logistic model for spouse diabetes concordance.
# ============================================================================

# Variable meaning in current pipeline:
# 1) diabetes: diabetes status indicator at a given wave (not incidence)
# 2) dm_incident: first incident diabetes event among baseline diabetes-free only

# SAS-equivalent analysis uses pair-level baseline diabetes status in wide format.
baseline_couples <- analytic_df_wide %>%
  dplyr::filter(fup == 0) %>%
  dplyr::select(hhid, carrs, female_dm_biomarker, male_dm_biomarker) %>%
  dplyr::filter(!is.na(female_dm_biomarker), !is.na(male_dm_biomarker))

# -----------------------------------------------------------------------------
# 1) Frequency table: female_dm_biomarker * male_dm_biomarker
# -----------------------------------------------------------------------------
dm_xtab <- with(baseline_couples, table(female_dm_biomarker, male_dm_biomarker))
print(dm_xtab)

# Relative risk / odds ratio analogue for 2x2 table
# fisher.test returns an exact odds ratio and confidence interval.
dm_fisher <- fisher.test(dm_xtab)
print(dm_fisher)

# -----------------------------------------------------------------------------
# 2) Logistic model: male_dm_biomarker (1) ~ female_dm_biomarker
# -----------------------------------------------------------------------------
dm_logit <- glm(male_dm_biomarker ~ female_dm_biomarker, data = baseline_couples, family = binomial())

coef_est <- coef(dm_logit)["female_dm_biomarker"]
coef_se <- sqrt(vcov(dm_logit)["female_dm_biomarker", "female_dm_biomarker"])
coef_p <- summary(dm_logit)$coefficients["female_dm_biomarker", "Pr(>|z|)"]

logit_or_tbl <- tibble::tibble(
  term = "female_dm_biomarker (1 vs 0)",
  estimate_log_odds = coef_est,
  odds_ratio = exp(coef_est),
  conf_low_95 = exp(coef_est - 1.96 * coef_se),
  conf_high_95 = exp(coef_est + 1.96 * coef_se),
  p_value = coef_p
)

print(summary(dm_logit))
print(logit_or_tbl)

# ============================================================================
# 3) dmcat2 concordance analysis (dmcat2 equals diabetes indicator)
# ============================================================================

dyads_dmcat2 <- baseline_couples

# Cross-tab: female_dm_biomarker * male_dm_biomarker
dmcat2_xtab <- with(dyads_dmcat2, table(female_dm_biomarker, male_dm_biomarker))
print(dmcat2_xtab)

# Relative risk / odds ratio analogue
dmcat2_fisher <- fisher.test(dmcat2_xtab)
print(dmcat2_fisher)

# Agreement metrics and McNemar test
agree_n <- sum(diag(dmcat2_xtab))
agree_total <- sum(dmcat2_xtab)
overall_agreement <- agree_n / agree_total
print(overall_agreement)

dmcat2_mcnemar <- mcnemar.test(dmcat2_xtab)
print(dmcat2_mcnemar)

# Chi-square association test
dmcat2_chisq <- chisq.test(dmcat2_xtab, correct = FALSE)
print(dmcat2_chisq)

# Binary logistic: male_dm_biomarker (1) ~ female_dm_biomarker
dmcat2_logit <- glm(male_dm_biomarker ~ female_dm_biomarker, data = dyads_dmcat2, family = binomial())

coef_est2 <- coef(dmcat2_logit)["female_dm_biomarker"]
coef_se2 <- sqrt(vcov(dmcat2_logit)["female_dm_biomarker", "female_dm_biomarker"])
coef_p2 <- summary(dmcat2_logit)$coefficients["female_dm_biomarker", "Pr(>|z|)"]

logit_dmcat2_or_tbl <- tibble::tibble(
  term = "female_dm_biomarker (1 vs 0)",
  estimate_log_odds = coef_est2,
  odds_ratio = exp(coef_est2),
  conf_low_95 = exp(coef_est2 - 1.96 * coef_se2),
  conf_high_95 = exp(coef_est2 + 1.96 * coef_se2),
  p_value = coef_p2
)

print(summary(dmcat2_logit))
print(logit_dmcat2_or_tbl)

# -----------------------------------------------------------------------------
# 4) Table and heat maps for baseline diabetes concordance
# -----------------------------------------------------------------------------

dm_table_labeled <- with(
  baseline_couples,
  table(
    women = factor(female_dm_biomarker, levels = c(0, 1), labels = c("No diabetes", "Diabetes")),
    men = factor(male_dm_biomarker, levels = c(0, 1), labels = c("No diabetes", "Diabetes"))
  )
)

dm_row_pct <- 100 * prop.table(dm_table_labeled, margin = 1)
dm_col_pct <- 100 * prop.table(dm_table_labeled, margin = 2)

dm_cell_df <- as.data.frame(as.table(dm_table_labeled), stringsAsFactors = FALSE) %>%
  dplyr::rename(Women_Status = women, Men_Status = men, Count = Freq) %>%
  dplyr::mutate(
    overall_pct = 100 * Count / sum(Count),
    row_pct = dm_row_pct[cbind(as.character(Women_Status), as.character(Men_Status))],
    col_pct = dm_col_pct[cbind(as.character(Women_Status), as.character(Men_Status))],
    cell_text = sprintf(
      "%d (%.1f%%)\n(%.1f%%)\n(%.1f%%)",
      Count,
      overall_pct,
      row_pct,
      col_pct
    )
  )

dm_table_display <- dm_cell_df %>%
  dplyr::select(Women_Status, Men_Status, cell_text) %>%
  tidyr::pivot_wider(names_from = Men_Status, values_from = cell_text)

print(dm_table_display)
utils::write.csv(
  dm_cell_df,
  file = "cca/analysis/psdcan02_dm concordance table cells.csv",
  row.names = FALSE
)

