# ============================================================
# 16_high_specificity_persistence_and_boundary_fix.R
#
# Finasteride / Dutasteride and Incident Meningioma
#
# FINAL targeted analysis patch:
#   1) Resolve the landmark-boundary bookkeeping issue.
#      Tumors occurring ON the landmark are treated as prevalent
#      at risk-set entry and excluded; incident outcomes must occur
#      STRICTLY AFTER the landmark.
#
#   2) Re-run the two key persistence findings using a
#      high-specificity intracranial meningioma endpoint:
#        - ICD-10-CM D32.0
#        - historical ICD-9-CM 225.2
#
#      Main interrogations:
#        a) strict sustained 1-y FINASTERIDE vs uro-alpha
#        b) persistent 3-y FINASTERIDE vs uro-alpha
#
#      All-5ARI versions are also reported for context.
#
#   3) Re-run the boundary-fixed primary and target-trial 1-y
#      models to verify exact cohort reconciliation.
#
#   4) Produce a conservative expanded BH-FDR table containing
#      the Stage-12 exploratory family plus these new
#      high-specificity persistence analyses.
#
# NO RAW CSV SCANS.
# Uses:
#   method_dx_events                         (Stage 14)
#   analysis_cohort_primary_revised          (Stage 05/06)
#   analysis_cohort_tte_index_itt_lag1       (Stage 14)
#   analysis_cohort_sustained_1y_strict_util (Stage 11)
#   analysis_cohort_persistent_3y_util        (Stage 11)
# ============================================================

run_stage_16 <- function() {

  options(stringsAsFactors = FALSE, scipen = 999)

  required_pkgs <- c(
    "DBI",
    "duckdb",
    "data.table",
    "survival",
    "splines"
  )

  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0L) {
    install.packages(missing_pkgs)
  }

  suppressPackageStartupMessages({
    library(DBI)
    library(duckdb)
    library(data.table)
    library(survival)
    library(splines)
  })

  # ----------------------------------------------------------
  # 0. PATHS
  # ----------------------------------------------------------
  R_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"
  WORK_DIR <- "C:/TriNetX/bph_m_work"
  DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")

  OUT_DIR <- file.path(
    WORK_DIR,
    "analysis_stage16_high_specificity"
  )

  QC_DIR <- file.path(
    WORK_DIR,
    "qc_stage16_high_specificity"
  )

  SIGNAL_DIR <- file.path(
    WORK_DIR,
    "analysis_signal_interrogation"
  )

  METHOD_DIR <- file.path(
    WORK_DIR,
    "analysis_methodologic_validity"
  )

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

  source(
    file.path(R_DIR, "final_analysis_helpers.R")
  )

  N_THREADS <- min(
    24L,
    max(4L, parallel::detectCores(logical = TRUE) - 6L)
  )
  MEMORY_LIMIT <- "48GB"

  con <- dbConnect(
    duckdb::duckdb(),
    dbdir = DB_PATH,
    read_only = FALSE
  )

  on.exit({
    try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }, add = TRUE)

  dbExecute(con, sprintf("PRAGMA threads=%d;", N_THREADS))
  dbExecute(con, sprintf("PRAGMA memory_limit='%s';", MEMORY_LIMIT))

  needed <- c(
    "method_dx_events",
    "analysis_cohort_primary_revised",
    "analysis_cohort_tte_index_itt_lag1",
    "analysis_cohort_sustained_1y_strict_util",
    "analysis_cohort_persistent_3y_util"
  )

  missing <- setdiff(needed, dbListTables(con))

  if (length(missing) > 0L) {
    stop(
      "Missing prerequisite tables. Run Stages 11-15 first. Missing: ",
      paste(missing, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 1. BUILD LANDMARK-BOUNDARY-FIXED TABLES
  # ----------------------------------------------------------
  #
  # Risk-set convention:
  #   - any broad meningeal tumor on/before landmark => excluded
  #   - outcome must occur STRICTLY AFTER landmark
  #
  # Existing source cohorts already applied the earlier eligibility
  # criteria; therefore this patch should ordinarily remove only
  # exact-boundary cases.
  #
  build_boundary_fixed <- function(source_table, output_table) {

    sql <- paste0(
"
CREATE OR REPLACE TABLE ", output_table, " AS

WITH outcome_dates AS (

  SELECT
    a.patient_id,

    MIN(m.dx_date) FILTER (
      WHERE m.is_broad_meningeal=1
        AND m.dx_date <= a.landmark_date
    ) AS broad_meningeal_on_or_before_landmark,

    MIN(m.dx_date) FILTER (
      WHERE m.is_strict_meningioma=1
        AND m.dx_date > a.landmark_date
        AND m.dx_date <= a.last_encounter_date
    ) AS strict_event_after_landmark,

    MIN(m.dx_date) FILTER (
      WHERE m.is_specific_intracranial=1
        AND m.dx_date > a.landmark_date
        AND m.dx_date <= a.last_encounter_date
    ) AS specific_event_after_landmark

  FROM ", source_table, " a

  LEFT JOIN method_dx_events m
    ON m.patient_id=a.patient_id

  GROUP BY a.patient_id
)

SELECT
  a.*,

  o.broad_meningeal_on_or_before_landmark,

  CASE
    WHEN o.strict_event_after_landmark IS NOT NULL
    THEN 1 ELSE 0
  END AS meningioma_event_strict_boundary,

  o.strict_event_after_landmark
    AS incident_meningioma_date_strict_boundary,

  DATE_DIFF(
    'day',
    a.landmark_date,
    CASE
      WHEN o.strict_event_after_landmark IS NOT NULL
      THEN o.strict_event_after_landmark
      ELSE a.last_encounter_date
    END
  ) AS followup_days_strict_boundary,

  CASE
    WHEN o.specific_event_after_landmark IS NOT NULL
    THEN 1 ELSE 0
  END AS meningioma_event_specific_boundary,

  o.specific_event_after_landmark
    AS incident_meningioma_date_specific_boundary,

  DATE_DIFF(
    'day',
    a.landmark_date,
    CASE
      WHEN o.specific_event_after_landmark IS NOT NULL
      THEN o.specific_event_after_landmark
      ELSE a.last_encounter_date
    END
  ) AS followup_days_specific_boundary

FROM ", source_table, " a

JOIN outcome_dates o USING (patient_id)

WHERE o.broad_meningeal_on_or_before_landmark IS NULL
  AND a.last_encounter_date > a.landmark_date;
"
    )

    dbExecute(con, sql)
  }

  build_boundary_fixed(
    "analysis_cohort_primary_revised",
    "analysis_cohort_primary_boundary_fixed"
  )

  build_boundary_fixed(
    "analysis_cohort_tte_index_itt_lag1",
    "analysis_cohort_tte_lag1_boundary_fixed"
  )

  build_boundary_fixed(
    "analysis_cohort_sustained_1y_strict_util",
    "analysis_cohort_sustained_1y_strict_boundary_fixed"
  )

  build_boundary_fixed(
    "analysis_cohort_persistent_3y_util",
    "analysis_cohort_persistent_3y_boundary_fixed"
  )

  # ----------------------------------------------------------
  # 2. BOUNDARY QC
  # ----------------------------------------------------------
  cohort_map <- data.table(
    original = c(
      "analysis_cohort_primary_revised",
      "analysis_cohort_tte_index_itt_lag1",
      "analysis_cohort_sustained_1y_strict_util",
      "analysis_cohort_persistent_3y_util"
    ),
    fixed = c(
      "analysis_cohort_primary_boundary_fixed",
      "analysis_cohort_tte_lag1_boundary_fixed",
      "analysis_cohort_sustained_1y_strict_boundary_fixed",
      "analysis_cohort_persistent_3y_boundary_fixed"
    ),
    label = c(
      "Primary 1-y landmark",
      "Target-trial day-0 assignment, 1-y landmark",
      "Strict sustained 1-y",
      "Persistent 3-y"
    )
  )

  boundary_qc <- rbindlist(
    lapply(
      seq_len(nrow(cohort_map)),
      function(i) {

        old_tbl <- cohort_map$original[i]
        new_tbl <- cohort_map$fixed[i]
        lab <- cohort_map$label[i]

        old <- as.data.table(
          dbGetQuery(
            con,
            paste0(
              "
              SELECT
                COUNT(*) AS n_original,
                SUM(meningioma_event_strict) AS strict_events_original
              FROM ",
              old_tbl,
              ";"
            )
          )
        )

        new <- as.data.table(
          dbGetQuery(
            con,
            paste0(
              "
              SELECT
                COUNT(*) AS n_boundary_fixed,
                SUM(meningioma_event_strict_boundary)
                  AS strict_events_boundary_fixed,
                SUM(meningioma_event_specific_boundary)
                  AS specific_events_boundary_fixed
              FROM ",
              new_tbl,
              ";"
            )
          )
        )

        data.table(
          cohort = lab,
          n_original = old$n_original,
          n_boundary_fixed = new$n_boundary_fixed,
          patients_removed_at_landmark =
            old$n_original - new$n_boundary_fixed,
          strict_events_original =
            old$strict_events_original,
          strict_events_boundary_fixed =
            new$strict_events_boundary_fixed,
          specific_events_boundary_fixed =
            new$specific_events_boundary_fixed
        )
      }
    ),
    fill = TRUE
  )

  fwrite(
    boundary_qc,
    file.path(
      QC_DIR,
      "16_landmark_boundary_reconciliation.csv"
    )
  )

  # Exact boundary patients, for audit only.
  boundary_audit <- rbindlist(
    lapply(
      seq_len(nrow(cohort_map)),
      function(i) {

        old_tbl <- cohort_map$original[i]
        lab <- cohort_map$label[i]

        as.data.table(
          dbGetQuery(
            con,
            paste0(
"
SELECT DISTINCT
  '", gsub("'", "''", lab), "' AS cohort,
  a.patient_id,
  a.exposure,
  a.index_date,
  a.landmark_date,
  m.dx_date AS boundary_tumor_date,
  m.code_system,
  m.code_norm
FROM ", old_tbl, " a
JOIN method_dx_events m
  ON m.patient_id=a.patient_id
WHERE m.is_broad_meningeal=1
  AND m.dx_date = a.landmark_date
ORDER BY a.patient_id;
"
            )
          )
        )
      }
    ),
    fill = TRUE
  )

  fwrite(
    boundary_audit,
    file.path(
      QC_DIR,
      "16_exact_landmark_tumor_audit.csv"
    )
  )

  # ----------------------------------------------------------
  # 3. READ/MODEL HELPERS
  # ----------------------------------------------------------
  common_cols <- c(
    "patient_id",
    "exposure",
    "index_drug",
    "index_finasteride_strength",
    "index_dutasteride_strength",
    "index_date",
    "age_at_index",
    "race",
    "ethnicity",
    "marital_status",
    "patient_regional_location",
    "patient_source_id",
    "hypertension",
    "diabetes",
    "obesity",
    "ckd",
    "smoking",
    "cad",
    "heart_failure",
    "cerebrovascular_disease",
    "hyperlipidemia",
    "osteoporosis",
    "alopecia",
    "urinary_retention",
    "uti",
    "hematuria",
    "bladder_stone",
    "n_encounters_prior_365",
    "n_encounters_prior_1095",
    "n_ambulatory_prior_365",
    "n_ed_prior_365",
    "n_inpatient_prior_365"
  )

  util_cols_1y <- c(
    "lm_encounters_total",
    "lm_encounter_days",
    "lm_ambulatory",
    "lm_ed",
    "lm_inpatient",
    "lm_y1_encounter_days"
  )

  util_cols_3y <- c(
    util_cols_1y,
    "lm_y2_encounter_days",
    "lm_y3_encounter_days"
  )

  read_table_cols <- function(tbl, extra_cols) {

    fields <- dbListFields(con, tbl)

    use <- unique(
      c(common_cols, extra_cols)
    )

    use <- use[use %in% fields]

    as.data.table(
      dbGetQuery(
        con,
        paste0(
          "SELECT ",
          paste(use, collapse = ", "),
          " FROM ",
          tbl,
          ";"
        )
      )
    )
  }

  add_landmark_logs <- function(d, three_year = FALSE) {

    d <- copy(d)

    vars <- c(
      "lm_encounters_total",
      "lm_encounter_days",
      "lm_ambulatory",
      "lm_ed",
      "lm_inpatient",
      "lm_y1_encounter_days"
    )

    if (three_year) {
      vars <- c(
        vars,
        "lm_y2_encounter_days",
        "lm_y3_encounter_days"
      )
    }

    for (v in vars) {
      if (v %in% names(d)) {
        d[
          ,
          paste0("log_", v) :=
            log1p(as.numeric(get(v)))
        ]
      }
    }

    d
  }

  lm_numeric_1y <- c(
    "log_lm_encounters_total",
    "log_lm_encounter_days",
    "log_lm_ambulatory",
    "log_lm_ed",
    "log_lm_inpatient"
  )

  lm_numeric_3y <- c(
    "log_lm_encounters_total",
    "log_lm_ambulatory",
    "log_lm_ed",
    "log_lm_inpatient",
    "log_lm_y1_encounter_days",
    "log_lm_y2_encounter_days",
    "log_lm_y3_encounter_days"
  )

  model_rows <- list()
  balance_rows <- list()
  incidence_rows <- list()
  diag_rows <- list()

  store <- function(x) {
    model_rows[[length(model_rows) + 1L]] <<- x$result
    balance_rows[[length(balance_rows) + 1L]] <<- x$balance
    incidence_rows[[length(incidence_rows) + 1L]] <<- x$incidence
    diag_rows[[length(diag_rows) + 1L]] <<- x$diagnostics
  }

  # ----------------------------------------------------------
  # 4. RE-RUN BOUNDARY-FIXED ANCHOR MODELS
  # ----------------------------------------------------------
  primary_fixed <- read_table_cols(
    "analysis_cohort_primary_boundary_fixed",
    c(
      "meningioma_event_strict_boundary",
      "followup_days_strict_boundary"
    )
  )

  primary_fixed[
    ,
    trt := as.integer(exposure == "5ARI")
  ]

  x <- run_ato_generic(
    primary_fixed,
    "Boundary-fixed primary 1-y landmark",
    "meningioma_event_strict_boundary",
    "followup_days_strict_boundary"
  )

  store(x)

  tte_fixed <- read_table_cols(
    "analysis_cohort_tte_lag1_boundary_fixed",
    c(
      "meningioma_event_strict_boundary",
      "followup_days_strict_boundary"
    )
  )

  tte_fixed[
    ,
    trt := as.integer(exposure == "5ARI")
  ]

  x <- run_ato_generic(
    tte_fixed,
    "Boundary-fixed target-trial day-0 assignment, 1-y landmark",
    "meningioma_event_strict_boundary",
    "followup_days_strict_boundary"
  )

  store(x)

  # ----------------------------------------------------------
  # 5. STRICT SUSTAINED 1-Y:
  #    HIGH-SPECIFICITY OUTCOME
  # ----------------------------------------------------------
  s1s <- read_table_cols(
    "analysis_cohort_sustained_1y_strict_boundary_fixed",
    c(
      "meningioma_event_strict_boundary",
      "followup_days_strict_boundary",
      "meningioma_event_specific_boundary",
      "followup_days_specific_boundary",
      util_cols_1y
    )
  )

  s1s <- add_landmark_logs(s1s)

  # All 5ARI vs uro-alpha, context only.
  s1s_all <- copy(s1s)
  s1s_all[
    ,
    trt := as.integer(exposure == "5ARI")
  ]

  x <- run_ato_generic(
    s1s_all,
    "High-specificity: strict sustained 1-y, all 5ARI vs uro-alpha",
    "meningioma_event_specific_boundary",
    "followup_days_specific_boundary",
    extra_numeric = lm_numeric_1y
  )

  store(x)

  # Finasteride-specific key interrogation.
  s1s_fin <- s1s[
    exposure == "URO_ALPHA" |
    (
      exposure == "5ARI" &
      index_drug == "finasteride"
    )
  ]

  s1s_fin[
    ,
    trt := as.integer(exposure == "5ARI")
  ]

  x <- run_ato_generic(
    s1s_fin,
    "High-specificity: strict sustained 1-y, finasteride vs uro-alpha",
    "meningioma_event_specific_boundary",
    "followup_days_specific_boundary",
    extra_numeric = lm_numeric_1y
  )

  store(x)

  # Also rerun the original strict endpoint after exact boundary repair,
  # so any change in the earlier positive estimate can be quantified.
  x <- run_ato_generic(
    s1s_fin,
    "Boundary-fixed strict endpoint: strict sustained 1-y, finasteride vs uro-alpha",
    "meningioma_event_strict_boundary",
    "followup_days_strict_boundary",
    extra_numeric = lm_numeric_1y
  )

  store(x)

  # ----------------------------------------------------------
  # 6. PERSISTENT 3-Y:
  #    HIGH-SPECIFICITY OUTCOME
  # ----------------------------------------------------------
  p3 <- read_table_cols(
    "analysis_cohort_persistent_3y_boundary_fixed",
    c(
      "meningioma_event_strict_boundary",
      "followup_days_strict_boundary",
      "meningioma_event_specific_boundary",
      "followup_days_specific_boundary",
      util_cols_3y
    )
  )

  p3 <- add_landmark_logs(
    p3,
    three_year = TRUE
  )

  # All 5ARI vs uro-alpha, context only.
  p3_all <- copy(p3)
  p3_all[
    ,
    trt := as.integer(exposure == "5ARI")
  ]

  x <- run_ato_generic(
    p3_all,
    "High-specificity: persistent 3-y, all 5ARI vs uro-alpha",
    "meningioma_event_specific_boundary",
    "followup_days_specific_boundary",
    extra_numeric = lm_numeric_3y
  )

  store(x)

  # Finasteride-specific key interrogation.
  p3_fin <- p3[
    exposure == "URO_ALPHA" |
    (
      exposure == "5ARI" &
      index_drug == "finasteride"
    )
  ]

  p3_fin[
    ,
    trt := as.integer(exposure == "5ARI")
  ]

  x <- run_ato_generic(
    p3_fin,
    "High-specificity: persistent 3-y, finasteride vs uro-alpha",
    "meningioma_event_specific_boundary",
    "followup_days_specific_boundary",
    extra_numeric = lm_numeric_3y
  )

  store(x)

  # Boundary-fixed original endpoint.
  x <- run_ato_generic(
    p3_fin,
    "Boundary-fixed strict endpoint: persistent 3-y, finasteride vs uro-alpha",
    "meningioma_event_strict_boundary",
    "followup_days_strict_boundary",
    extra_numeric = lm_numeric_3y
  )

  store(x)

  # ----------------------------------------------------------
  # 7. SAVE MODEL OUTPUTS
  # ----------------------------------------------------------
  models <- rbindlist(
    model_rows,
    fill = TRUE
  )

  balances <- rbindlist(
    balance_rows,
    fill = TRUE
  )

  incidence <- rbindlist(
    incidence_rows,
    fill = TRUE
  )

  diagnostics <- rbindlist(
    diag_rows,
    fill = TRUE
  )

  fwrite(
    models,
    file.path(
      OUT_DIR,
      "16_model_results.csv"
    )
  )

  fwrite(
    balances,
    file.path(
      OUT_DIR,
      "16_balance_all.csv"
    )
  )

  fwrite(
    balance_summary(balances),
    file.path(
      OUT_DIR,
      "16_balance_summary.csv"
    )
  )

  fwrite(
    incidence,
    file.path(
      OUT_DIR,
      "16_incidence_rates.csv"
    )
  )

  fwrite(
    diagnostics,
    file.path(
      OUT_DIR,
      "16_predictor_diagnostics.csv"
    )
  )

  # ----------------------------------------------------------
  # 8. CONSERVATIVE EXPANDED EXPLORATORY FDR
  # ----------------------------------------------------------
  #
  # Stage 12 already had a 12-model exploratory family.
  # Here we add the FOUR high-specificity persistence models.
  # The two boundary-only reruns and anchor models are not new
  # hypotheses and are therefore not added to this family.
  #
  stage12_file <- file.path(
    SIGNAL_DIR,
    "12_signal_model_results.csv"
  )

  if (file.exists(stage12_file)) {

    old <- fread(stage12_file)

    old_fdr <- old[
      ,
      .(
        analysis,
        weighted_HR,
        weighted_CI_low,
        weighted_CI_high,
        weighted_p,
        family_source = "Stage 12 exploratory"
      )
    ]

    new_fdr <- models[
      grepl(
        "^High-specificity:",
        analysis
      ),
      .(
        analysis,
        weighted_HR,
        weighted_CI_low,
        weighted_CI_high,
        weighted_p,
        family_source = "Stage 16 high-specificity sensitivity"
      )
    ]

    expanded <- rbindlist(
      list(
        old_fdr,
        new_fdr
      ),
      fill = TRUE
    )

    expanded[
      ,
      q_BH_expanded := p.adjust(
        weighted_p,
        method = "BH"
      )
    ]

    expanded[
      ,
      `:=`(
        raw_p_lt_0_05 =
          as.integer(weighted_p < 0.05),
        q_BH_expanded_lt_0_05 =
          as.integer(q_BH_expanded < 0.05)
      )
    ]

    fwrite(
      expanded,
      file.path(
        OUT_DIR,
        "16_expanded_exploratory_FDR.csv"
      )
    )
  }

  # ----------------------------------------------------------
  # 9. COMPARE WITH STAGE-15 IPCW COHORT COUNT
  # ----------------------------------------------------------
  ipcw_file <- file.path(
    METHOD_DIR,
    "15_ipcw_results.csv"
  )

  if (file.exists(ipcw_file)) {

    ipcw <- fread(ipcw_file)

    fixed_primary_counts <- models[
      analysis == "Boundary-fixed primary 1-y landmark",
      .(
        boundary_fixed_primary_n = n,
        boundary_fixed_primary_events = events
      )
    ]

    ipcw_compare <- cbind(
      fixed_primary_counts,
      ipcw[
        grepl(
          "1st/99th",
          analysis
        ),
        .(
          ipcw_n = n_patients,
          ipcw_events = events,
          ipcw_HR = HR,
          ipcw_CI_low = CI_low,
          ipcw_CI_high = CI_high,
          ipcw_p = p
        )
      ]
    )

    fwrite(
      ipcw_compare,
      file.path(
        QC_DIR,
        "16_boundary_fixed_vs_ipcw_reconciliation.csv"
      )
    )
  }

  # ----------------------------------------------------------
  # 10. COMPACT REVIEW FILE
  # ----------------------------------------------------------
  review_labels <- c(
    "Boundary-fixed primary 1-y landmark",
    "Boundary-fixed target-trial day-0 assignment, 1-y landmark",
    "Boundary-fixed strict endpoint: strict sustained 1-y, finasteride vs uro-alpha",
    "High-specificity: strict sustained 1-y, finasteride vs uro-alpha",
    "Boundary-fixed strict endpoint: persistent 3-y, finasteride vs uro-alpha",
    "High-specificity: persistent 3-y, finasteride vs uro-alpha"
  )

  review <- models[
    analysis %in% review_labels,
    .(
      analysis,
      n,
      n_treated,
      n_control,
      events,
      events_treated,
      events_control,
      weighted_HR,
      weighted_CI_low,
      weighted_CI_high,
      weighted_p,
      ESS_treated,
      ESS_control
    )
  ]

  review[
    ,
    order_id := match(
      analysis,
      review_labels
    )
  ]

  setorder(review, order_id)
  review[, order_id := NULL]

  fwrite(
    review,
    file.path(
      OUT_DIR,
      "16_REVIEW_FIRST.csv"
    )
  )

  cat("\n================ BOUNDARY RECONCILIATION ============\n")
  print(boundary_qc)

  cat("\n================ STAGE 16 KEY RESULTS ================\n")
  print(review)

  cat("\nStage 16 complete.\n")
  cat("Analysis directory:", OUT_DIR, "\n")
  cat("QC directory:", QC_DIR, "\n")

  cat(
    "\nFIRST FILE TO SEND BACK:\n",
    file.path(
      OUT_DIR,
      "16_REVIEW_FIRST.csv"
    ),
    "\n",
    sep = ""
  )
}

run_stage_16()
