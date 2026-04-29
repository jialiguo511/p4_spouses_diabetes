rm(list=ls());gc();source(".Rprofile")

source("cca/preprocessing/archive/psdcpre03_dm incidence df.R")

analytic_df_wide <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/archive/psdcpre04_analytic dataset wide.RDS"))

# ============================================================================
# Table 3: Risk of diabetes over follow-up by partner baseline diabetes status
# among individuals without diabetes at baseline.
# ============================================================================

# One row per couple at baseline for partner baseline status.
baseline_dyads <- analytic_df_wide %>%
  dplyr::filter(fup == 0) %>%
  dplyr::select(
    hhid,
    carrs,
    female_pid,
    male_pid,
    female_baseline_dm,
    male_baseline_dm
  ) %>%
  dplyr::distinct()

# Person-level follow-up diabetes outcome among baseline diabetes-free participants.
followup_outcome <- dm_event_or_censor_df %>%
  dplyr::select(carrs, pid, dm_incident) %>%
  dplyr::distinct(carrs, pid, .keep_all = TRUE)

build_table_section <- function(data, baseline_label) {
  data %>%
    dplyr::mutate(
      partner_status = factor(
        partner_status,
        levels = c("No diabetes", "Diabetes")
      ),
      outcome_status = factor(
        outcome_status,
        levels = c("No diabetes", "Diabetes")
      )
    ) %>%
    dplyr::count(partner_status, outcome_status, name = "n") %>%
    dplyr::group_by(partner_status) %>%
    dplyr::mutate(
      partner_n = sum(n),
      pct = 100 * n / partner_n
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      baseline_group = baseline_label,
      partner_group = paste0(as.character(partner_status), " (n=", partner_n, ")")
    ) %>%
    dplyr::select(baseline_group, partner_group, outcome_status, pct) %>%
    tidyr::pivot_wider(
      names_from = outcome_status,
      values_from = pct
    ) %>%
    dplyr::rename(
      no_diabetes_pct = `No diabetes`,
      diabetes_pct = Diabetes
    ) %>%
    dplyr::mutate(
      no_diabetes_pct = sprintf("%.1f", no_diabetes_pct),
      diabetes_pct = sprintf("%.1f", diabetes_pct)
    )
}

# Women without diabetes at baseline: follow-up diabetes by husband's baseline status.
women_df <- baseline_dyads %>%
  dplyr::filter(female_baseline_dm == 0) %>%
  dplyr::left_join(
    followup_outcome %>%
      dplyr::rename(female_pid = pid, female_followup_dm = dm_incident),
    by = c("carrs", "female_pid")
  ) %>%
  dplyr::mutate(
    partner_status = dplyr::case_when(
      male_baseline_dm == 0 ~ "No diabetes",
      male_baseline_dm == 1 ~ "Diabetes",
      TRUE ~ NA_character_
    ),
    outcome_status = dplyr::case_when(
      female_followup_dm == 0 ~ "No diabetes",
      female_followup_dm == 1 ~ "Diabetes",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(partner_status), !is.na(outcome_status))

women_label <- paste0(
  "Women without diabetes at baseline (n=",
  nrow(women_df),
  ")"
)

women_section <- build_table_section(women_df, women_label)

# Men without diabetes at baseline: follow-up diabetes by wife's baseline status.
men_df <- baseline_dyads %>%
  dplyr::filter(male_baseline_dm == 0) %>%
  dplyr::left_join(
    followup_outcome %>%
      dplyr::rename(male_pid = pid, male_followup_dm = dm_incident),
    by = c("carrs", "male_pid")
  ) %>%
  dplyr::mutate(
    partner_status = dplyr::case_when(
      female_baseline_dm == 0 ~ "No diabetes",
      female_baseline_dm == 1 ~ "Diabetes",
      TRUE ~ NA_character_
    ),
    outcome_status = dplyr::case_when(
      male_followup_dm == 0 ~ "No diabetes",
      male_followup_dm == 1 ~ "Diabetes",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(partner_status), !is.na(outcome_status))

men_label <- paste0(
  "Men without diabetes at baseline (n=",
  nrow(men_df),
  ")"
)

men_section <- build_table_section(men_df, men_label)

table3_df <- dplyr::bind_rows(women_section, men_section) %>%
  dplyr::rename(
    `Individual's status at baseline` = baseline_group,
    `Individuals with following partner's status at baseline` = partner_group,
    `No diabetes, %` = no_diabetes_pct,
    `Diabetes, %` = diabetes_pct
  )

print(table3_df)

utils::write.csv(
  table3_df,
  file = "cca/analysis/psdcan03_followup_risk_by_spouse_baseline.csv",
  row.names = FALSE
)
