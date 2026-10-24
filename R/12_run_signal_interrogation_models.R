# ============================================================
# 12_run_signal_interrogation_models.R
#
# Tests whether the two exploratory positive signals survive:
#   1) adjustment for healthcare utilization measured DURING the
#      exposure-assessment landmark period;
#   2) drug-specific sustained/persistent comparisons;
#   3) homogeneous contemporary US-only exact-dose
#      finasteride-vs-dutasteride restrictions.
#
# These are targeted sensitivity analyses, not a new primary analysis.
# ============================================================

run_stage_12 <- function() {

  options(stringsAsFactors = FALSE, scipen = 999)

  required_pkgs <- c(
    "DBI", "duckdb", "data.table", "survival", "splines"
  )
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

  library(DBI)
  library(duckdb)
  library(data.table)
  library(survival)
  library(splines)

  R_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"
  WORK_DIR <- "C:/TriNetX/bph_m_work"
  DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")
  OUT_DIR <- file.path(WORK_DIR, "analysis_signal_interrogation")
  PREV_RESULTS <- file.path(
    WORK_DIR,
    "analysis_final",
    "10_final_model_results.csv"
  )

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  source(file.path(R_DIR, "final_analysis_helpers.R"))

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
    "analysis_cohort_sustained_1y_util",
    "analysis_cohort_sustained_1y_strict_util",
    "analysis_cohort_persistent_3y_util",
    "analysis_cohort_fd_us_2012_exact",
    "analysis_cohort_fd_us_2012_exact_age50",
    "analysis_cohort_fd_us_icd10_exact"
  )
  missing <- setdiff(needed, dbListTables(con))
  if (length(missing) > 0L) {
    stop(
      "Run Stage 11 first. Missing: ",
      paste(missing, collapse = ", ")
    )
  }

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
    use <- unique(c(common_cols, extra_cols))
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
        d[, paste0("log_", v) := log1p(as.numeric(get(v)))]
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

  # For the 3-y model, year-specific encounter-day variables let us balance
  # ongoing contact across the exposure-assessment period without using a
  # huge list of redundant raw counts.
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
  diag_rows <- list()
  incidence_rows <- list()

  store_result <- function(x) {
    model_rows[[length(model_rows) + 1L]] <<- x$result
    balance_rows[[length(balance_rows) + 1L]] <<- x$balance
    diag_rows[[length(diag_rows) + 1L]] <<- x$diagnostics
    incidence_rows[[length(incidence_rows) + 1L]] <<- x$incidence
  }

  safe_run <- function(
      d,
      label,
      extra_numeric = character(),
      min_treated = 300,
      min_control = 300,
      min_events = 20
  ) {

    if (length(unique(d$trt)) < 2L) {
      warning(label, ": only one treatment group; skipped.")
      return(NULL)
    }

    if (
      nrow(d[trt == 1]) < min_treated ||
      nrow(d[trt == 0]) < min_control ||
      sum(d$meningioma_event_strict, na.rm = TRUE) < min_events
    ) {
      warning(label, ": sparse sample/event count; skipped.")
      return(NULL)
    }

    run_ato_generic(
      d,
      label = label,
      event_col = "meningioma_event_strict",
      time_days_col = "followup_days_strict",
      extra_numeric = extra_numeric
    )
  }

  # ----------------------------------------------------------
  # 1. 1-Y SUSTAINED: ALL 5ARI, FINASTERIDE, DUTASTERIDE
  # ----------------------------------------------------------
  s1 <- read_table_cols(
    "analysis_cohort_sustained_1y_util",
    c(
      "meningioma_event_strict",
      "followup_days_strict",
      util_cols_1y
    )
  )
  s1 <- add_landmark_logs(s1)

  s1[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    s1,
    "Sustained 1-y + landmark utilization: all 5ARI vs uro-alpha",
    lm_numeric_1y
  )
  if (!is.null(x)) store_result(x)

  s1_fin <- s1[
    exposure == "URO_ALPHA" |
    (exposure == "5ARI" & index_drug == "finasteride")
  ]
  s1_fin[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    s1_fin,
    "Sustained 1-y + landmark utilization: finasteride vs uro-alpha",
    lm_numeric_1y
  )
  if (!is.null(x)) store_result(x)

  s1_dut <- s1[
    exposure == "URO_ALPHA" |
    (exposure == "5ARI" & index_drug == "dutasteride")
  ]
  s1_dut[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    s1_dut,
    "Sustained 1-y + landmark utilization: dutasteride vs uro-alpha",
    lm_numeric_1y,
    min_treated = 200
  )
  if (!is.null(x)) store_result(x)

  # ----------------------------------------------------------
  # 2. STRICT 1-Y SUSTAINED
  # ----------------------------------------------------------
  s1s <- read_table_cols(
    "analysis_cohort_sustained_1y_strict_util",
    c(
      "meningioma_event_strict",
      "followup_days_strict",
      util_cols_1y
    )
  )
  s1s <- add_landmark_logs(s1s)

  s1s[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    s1s,
    "Strict sustained 1-y + landmark utilization: all 5ARI vs uro-alpha",
    lm_numeric_1y
  )
  if (!is.null(x)) store_result(x)

  s1s_fin <- s1s[
    exposure == "URO_ALPHA" |
    (exposure == "5ARI" & index_drug == "finasteride")
  ]
  s1s_fin[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    s1s_fin,
    "Strict sustained 1-y + landmark utilization: finasteride vs uro-alpha",
    lm_numeric_1y
  )
  if (!is.null(x)) store_result(x)

  s1s_dut <- s1s[
    exposure == "URO_ALPHA" |
    (exposure == "5ARI" & index_drug == "dutasteride")
  ]
  s1s_dut[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    s1s_dut,
    "Strict sustained 1-y + landmark utilization: dutasteride vs uro-alpha",
    lm_numeric_1y,
    min_treated = 100,
    min_events = 15
  )
  if (!is.null(x)) store_result(x)

  # ----------------------------------------------------------
  # 3. 3-Y PERSISTENT
  # ----------------------------------------------------------
  p3 <- read_table_cols(
    "analysis_cohort_persistent_3y_util",
    c(
      "meningioma_event_strict",
      "followup_days_strict",
      util_cols_3y
    )
  )
  p3 <- add_landmark_logs(p3, three_year = TRUE)

  p3[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    p3,
    "Persistent 3-y + landmark utilization: all 5ARI vs uro-alpha",
    lm_numeric_3y
  )
  if (!is.null(x)) store_result(x)

  p3_fin <- p3[
    exposure == "URO_ALPHA" |
    (exposure == "5ARI" & index_drug == "finasteride")
  ]
  p3_fin[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    p3_fin,
    "Persistent 3-y + landmark utilization: finasteride vs uro-alpha",
    lm_numeric_3y
  )
  if (!is.null(x)) store_result(x)

  p3_dut <- p3[
    exposure == "URO_ALPHA" |
    (exposure == "5ARI" & index_drug == "dutasteride")
  ]
  p3_dut[, trt := as.integer(exposure == "5ARI")]
  x <- safe_run(
    p3_dut,
    "Persistent 3-y + landmark utilization: dutasteride vs uro-alpha",
    lm_numeric_3y,
    min_treated = 100,
    min_events = 15
  )
  if (!is.null(x)) store_result(x)

  # ----------------------------------------------------------
  # 4. HOMOGENEOUS FINASTERIDE VS DUTASTERIDE
  # ----------------------------------------------------------
  run_fd <- function(tbl, label) {
    d <- read_table_cols(
      tbl,
      c(
        "meningioma_event_strict",
        "followup_days_strict"
      )
    )
    d[, trt := as.integer(index_drug == "finasteride")]

    safe_run(
      d,
      label,
      extra_numeric = character(),
      min_treated = 300,
      min_control = 300,
      min_events = 20
    )
  }

  x <- run_fd(
    "analysis_cohort_fd_us_2012_exact",
    "Head-to-head homogeneous: US, 2012+, exact BPH doses, finasteride vs dutasteride"
  )
  if (!is.null(x)) store_result(x)

  x <- run_fd(
    "analysis_cohort_fd_us_2012_exact_age50",
    "Head-to-head homogeneous: US, 2012+, age >=50, exact BPH doses"
  )
  if (!is.null(x)) store_result(x)

  x <- run_fd(
    "analysis_cohort_fd_us_icd10_exact",
    "Head-to-head homogeneous: US, ICD-10 era, exact BPH doses"
  )
  if (!is.null(x)) store_result(x)

  # ----------------------------------------------------------
  # 5. SAVE RESULTS
  # ----------------------------------------------------------
  models <- rbindlist(model_rows, fill = TRUE)
  balances <- rbindlist(balance_rows, fill = TRUE)
  diagnostics <- rbindlist(diag_rows, fill = TRUE)
  incidence <- rbindlist(incidence_rows, fill = TRUE)

  fwrite(
    models,
    file.path(OUT_DIR, "12_signal_model_results.csv")
  )
  fwrite(
    balances,
    file.path(OUT_DIR, "12_signal_balance_all.csv")
  )
  fwrite(
    balance_summary(balances),
    file.path(OUT_DIR, "12_signal_balance_summary.csv")
  )
  fwrite(
    diagnostics,
    file.path(OUT_DIR, "12_signal_predictor_diagnostics.csv")
  )
  fwrite(
    incidence,
    file.path(OUT_DIR, "12_signal_incidence_rates.csv")
  )

  # ----------------------------------------------------------
  # 6. SIDE-BY-SIDE WITH PRIOR POSITIVE MODELS
  # ----------------------------------------------------------
  if (file.exists(PREV_RESULTS)) {

    old <- fread(PREV_RESULTS)

    old_keep <- old[
      grepl(
        "Sustained treatment|Persistent treatment|Head-to-head",
        analysis,
        ignore.case = TRUE
      )
    ]

    old_keep[, source_stage := "Stage 10"]
    new_keep <- copy(models)
    new_keep[, source_stage := "Stage 12"]

    comparison <- rbindlist(
      list(old_keep, new_keep),
      fill = TRUE
    )

    fwrite(
      comparison,
      file.path(OUT_DIR, "12_old_vs_targeted_models.csv")
    )
  }

  cat("\n================ TARGETED SIGNAL MODELS =============\n")
  print(models)

  cat("\n================ BALANCE SUMMARY ====================\n")
  print(balance_summary(balances))

  cat("\nStage 12 complete.\n")
  cat("Outputs:", OUT_DIR, "\n")
}

run_stage_12()
