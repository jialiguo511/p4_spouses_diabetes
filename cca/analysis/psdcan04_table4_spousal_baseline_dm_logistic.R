rm(list = ls()); gc(); source(".Rprofile")

# ============================================================================
# Table 4: Risk of diabetes diagnosis by spousal baseline diabetes status
# among individuals without diabetes at baseline.
#
# Models are fit separately by sex:
# 1) Age adjusted
# 2) Model 1: age + education + family history of diabetes
# 3) Model 2:
#    - Women: Model 1 + morbidity category + BMI category
#    - Men:   Model 1 + morbidity category + BMI category + alcohol + smoking
# ============================================================================

analytic_df <- readRDS(
  paste0(path_spouses_diabetes_folder, "/working/cca/preprocessing/psdcpre04_analytic dataset long.RDS")
)

analytic_df_wide <- readRDS(
  paste0(path_spouses_diabetes_folder, "/working/cca/preprocessing/psdcpre05_analytic dataset wide.RDS")
)

# One row per participant for incident diabetes outcome.
followup_outcome <- analytic_df %>%
  dplyr::select(carrs, pid, dm_incident) %>%
  dplyr::filter(!is.na(dm_incident)) %>%
  dplyr::distinct(carrs, pid, .keep_all = TRUE)

# One row per dyad at baseline with sex-specific covariates.
baseline_dyads <- analytic_df_wide %>%
  dplyr::filter(fup == 0) %>%
  dplyr::select(
    hhid,
    carrs,
    female_pid,
    male_pid,
    female_baseline_dm,
    male_baseline_dm,
    female_age,
    male_age,
    female_edu_category,
    male_edu_category,
    female_famhx_dm,
    male_famhx_dm,
    female_morbidity_category,
    male_morbidity_category,
    female_bmi_category,
    male_bmi_category,
    male_alc_overall,
    male_smk_overall
  ) %>%
  dplyr::distinct()

extract_partner_or <- function(model_obj, term_name = "partner_dm_statusDiabetes") {
  coef_est <- stats::coef(model_obj)[term_name]
  coef_se <- sqrt(stats::vcov(model_obj)[term_name, term_name])

  tibble::tibble(
    or = exp(coef_est),
    low = exp(coef_est - 1.96 * coef_se),
    high = exp(coef_est + 1.96 * coef_se)
  )
}

fmt_or <- function(or, low, high) {
  sprintf("%.2f (%.2f - %.2f)", or, low, high)
}

fit_sex_logit_models <- function(df, sex_label, include_behavior_terms = FALSE) {
  df_model <- df %>%
    dplyr::filter(
      baseline_dm == 0,
      !is.na(partner_baseline_dm),
      !is.na(dm_incident)
    ) %>%
    dplyr::mutate(
      dm_incident = as.integer(dm_incident),
      partner_dm_status = factor(
        partner_baseline_dm,
        levels = c(0, 1),
        labels = c("No diabetes", "Diabetes")
      ),
      edu_category = as.factor(edu_category),
      famhx_dm = as.factor(famhx_dm),
      morbidity_category = as.factor(morbidity_category),
      bmi_category = as.factor(bmi_category),
      alc_overall = as.factor(alc_overall),
      smk_overall = as.factor(smk_overall)
    )

  model_age <- stats::glm(
    dm_incident ~ partner_dm_status + age,
    data = df_model,
    family = stats::binomial()
  )

  model_1 <- stats::glm(
    dm_incident ~ partner_dm_status + age + edu_category + famhx_dm,
    data = df_model,
    family = stats::binomial()
  )

  if (include_behavior_terms) {
    model_2 <- stats::glm(
      dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
        morbidity_category + bmi_category + alc_overall + smk_overall,
      data = df_model,
      family = stats::binomial()
    )
  } else {
    model_2 <- stats::glm(
      dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
        morbidity_category + bmi_category,
      data = df_model,
      family = stats::binomial()
    )
  }

  age_est <- extract_partner_or(model_age)
  m1_est <- extract_partner_or(model_1)
  m2_est <- extract_partner_or(model_2)

  partner_n <- df_model %>%
    dplyr::count(partner_dm_status, name = "n")

  no_dm_n <- partner_n$n[partner_n$partner_dm_status == "No diabetes"]
  dm_n <- partner_n$n[partner_n$partner_dm_status == "Diabetes"]

  baseline_label <- paste0(
    sex_label,
    " without diabetes at baseline (n=",
    nrow(df_model),
    ")"
  )

  table_rows <- tibble::tibble(
    `Individual's status at baseline` = c(baseline_label, baseline_label),
    `Individuals with following partner's status at baseline` = c(
      paste0("No diabetes (n=", no_dm_n, ")"),
      paste0("Diabetes (n=", dm_n, ")")
    ),
    `Age adjusted OR (95% CI)` = c("1.00", fmt_or(age_est$or, age_est$low, age_est$high)),
    `Model 1 adjusted OR (95% CI)` = c("1.00", fmt_or(m1_est$or, m1_est$low, m1_est$high)),
    `Model 2 adjusted OR (95% CI)` = c("1.00", fmt_or(m2_est$or, m2_est$low, m2_est$high))
  )

  plot_rows <- tibble::tibble(
    sex = sex_label,
    model = c("Age adjusted", "Model 1", "Model 2"),
    or = c(age_est$or, m1_est$or, m2_est$or),
    low = c(age_est$low, m1_est$low, m2_est$low),
    high = c(age_est$high, m1_est$high, m2_est$high)
  )

  list(table_rows = table_rows, plot_rows = plot_rows)
}

women_df <- baseline_dyads %>%
  dplyr::transmute(
    hhid,
    carrs,
    pid = female_pid,
    baseline_dm = female_baseline_dm,
    partner_baseline_dm = male_baseline_dm,
    age = female_age,
    edu_category = female_edu_category,
    famhx_dm = female_famhx_dm,
    morbidity_category = female_morbidity_category,
    bmi_category = female_bmi_category,
    alc_overall = NA,
    smk_overall = NA
  ) %>%
  dplyr::left_join(followup_outcome, by = c("carrs", "pid"))

men_df <- baseline_dyads %>%
  dplyr::transmute(
    hhid,
    carrs,
    pid = male_pid,
    baseline_dm = male_baseline_dm,
    partner_baseline_dm = female_baseline_dm,
    age = male_age,
    edu_category = male_edu_category,
    famhx_dm = male_famhx_dm,
    morbidity_category = male_morbidity_category,
    bmi_category = male_bmi_category,
    alc_overall = male_alc_overall,
    smk_overall = male_smk_overall
  ) %>%
  dplyr::left_join(followup_outcome, by = c("carrs", "pid"))

women_results <- fit_sex_logit_models(
  df = women_df,
  sex_label = "Women",
  include_behavior_terms = FALSE
)

men_results <- fit_sex_logit_models(
  df = men_df,
  sex_label = "Men",
  include_behavior_terms = TRUE
)

table4_df <- dplyr::bind_rows(
  women_results$table_rows,
  men_results$table_rows
)

partner_or_plot_df <- dplyr::bind_rows(
  women_results$plot_rows,
  men_results$plot_rows
) %>%
  dplyr::mutate(
    model = factor(model, levels = c("Age adjusted", "Model 1", "Model 2")),
    sex = factor(sex, levels = c("Women", "Men"))
  )

print(table4_df)
print(partner_or_plot_df)

utils::write.csv(
  table4_df,
  file = "cca/analysis/psdcan04_table4_spousal_baseline_dm_logistic.csv",
  row.names = FALSE
)

utils::write.csv(
  partner_or_plot_df,
  file = "cca/analysis/psdcan04_table4_spousal_baseline_dm_logistic_partner_or.csv",
  row.names = FALSE
)

table4_plot <- ggplot2::ggplot(
  partner_or_plot_df,
  ggplot2::aes(
    x = model,
    y = or,
    ymin = low,
    ymax = high,
    color = sex,
    group = sex
  )
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  ggplot2::geom_pointrange(
    position = ggplot2::position_dodge(width = 0.35),
    linewidth = 0.45
  ) +
  ggplot2::labs(
    title = "Risk of diabetes diagnosis by spousal baseline diabetes status",
    subtitle = "Odds ratios for partner diabetes (Diabetes vs No diabetes)",
    x = "Model",
    y = "Odds Ratio (95% CI)",
    color = "Individual"
  ) +
  ggplot2::theme_bw(base_size = 11)

print(table4_plot)

ggplot2::ggsave(
  filename = "cca/analysis/psdcan04_table4_spousal_baseline_dm_logistic_plot.png",
  plot = table4_plot,
  width = 8,
  height = 5,
  dpi = 300
)
