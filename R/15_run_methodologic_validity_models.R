# ============================================================
# 15_run_methodologic_validity_models.R
#
# Runs the targeted methodological validity analyses built in Stage 14:
#
#   1) Target-trial-aligned day-0 ITT-like analyses
#      - zero lag
#      - 1-y biologic-latency landmark
#
#   2) EHR-continuity IPCW sensitivity
#      - baseline treatment overlap weights
#      - pooled-logistic censoring model in 180-d intervals
#      - time-updated prior-interval healthcare utilization
#      - 1/99 and 2.5/97.5 percentile truncation
#
#   3) Prespecified negative-control outcome analyses
#
#   4) Modern-era class-level sensitivities
#      - 2012 onward
#      - ICD-10 era (2015-10-01 onward)
#
#   5) Benjamini-Hochberg FDR for the Stage-12 exploratory family
#
# Outputs are deliberately separate from the finalized primary analysis.
# ============================================================

run_stage_15 <- function() {

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

  # speedglm is optional but strongly preferred for the several-million-row
  # pooled-logistic censoring model.
  if (!requireNamespace("speedglm", quietly = TRUE)) {
    try(
      install.packages("speedglm"),
      silent = TRUE
    )
  }

  # ----------------------------------------------------------
  # 0. PATHS / SETTINGS
  # ----------------------------------------------------------
  R_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"
  WORK_DIR <- "C:/TriNetX/bph_m_work"
  DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")
  OUT_DIR <- file.path(
    WORK_DIR,
    "analysis_methodologic_validity"
  )
  QC_DIR <- file.path(
    WORK_DIR,
    "qc_methodologic_validity"
  )
  SIGNAL_DIR <- file.path(
    WORK_DIR,
    "analysis_signal_interrogation"
  )

  dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )

  source(
    file.path(R_DIR, "final_analysis_helpers.R")
  )

  N_THREADS <- min(
    24L,
    max(4L, parallel::detectCores(logical = TRUE) - 6L)
  )
  MEMORY_LIMIT <- "48GB"

  IPCW_INTERVAL_DAYS <- 180L
  ADMIN_GRACE_DAYS <- 180L

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

  required_tables <- c(
    "analysis_cohort_tte_index_itt_zero",
    "analysis_cohort_tte_index_itt_lag1",
    "analysis_cohort_nco_cranial_nerve",
    "analysis_cohort_nco_spinal_cord",
    "analysis_cohort_primary_revised",
    "method_primary_180d_util",
    "candidate_observation"
  )

  missing_tables <- setdiff(
    required_tables,
    dbListTables(con)
  )

  if (length(missing_tables) > 0L) {
    stop(
      "Run Stage 14 first. Missing: ",
      paste(missing_tables, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 1. SHARED TABLE-READING HELPERS
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

  model_rows <- list()
  balance_rows <- list()
  diag_rows <- list()
  incidence_rows <- list()
  ph_rows <- list()

  store_result <- function(x) {

    model_rows[[length(model_rows) + 1L]] <<- x$result
    balance_rows[[length(balance_rows) + 1L]] <<- x$balance
    diag_rows[[length(diag_rows) + 1L]] <<- x$diagnostics
    incidence_rows[[length(incidence_rows) + 1L]] <<- x$incidence

    if (!is.null(x$zph)) {
      ph_rows[[length(ph_rows) + 1L]] <<- x$zph
    }
  }

  # ----------------------------------------------------------
  # 2. TARGET-TRIAL-ALIGNED DAY-0 ANALYSES
  # ----------------------------------------------------------
  tte0 <- read_table_cols(
    "analysis_cohort_tte_index_itt_zero",
    c(
      "meningioma_event_strict",
      "followup_days_strict"
    )
  )

  tte0[, trt := as.integer(exposure == "5ARI")]

  res <- run_ato_generic(
    tte0,
    "Target-trial aligned: day-0 ITT-like initiation effect, zero lag",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  tte1 <- read_table_cols(
    "analysis_cohort_tte_index_itt_lag1",
    c(
      "meningioma_event_strict",
      "followup_days_strict"
    )
  )

  tte1[, trt := as.integer(exposure == "5ARI")]

  res <- run_ato_generic(
    tte1,
    "Target-trial aligned: day-0 assignment, 1-y biologic-latency landmark",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  # ----------------------------------------------------------
  # 3. MODERN-ERA CLASS-LEVEL SENSITIVITIES
  # ----------------------------------------------------------
  primary <- read_table_cols(
    "analysis_cohort_primary_revised",
    c(
      "meningioma_event_strict",
      "followup_days_strict",
      "last_encounter_date",
      "landmark_date",
      "incident_meningioma_date_strict"
    )
  )

  primary[, index_date := as.Date(index_date)]
  primary[, trt := as.integer(exposure == "5ARI")]

  modern_2012 <- primary[
    index_date >= as.Date("2012-01-01")
  ]

  res <- run_ato_generic(
    modern_2012,
    "Modern-era primary cohort: index 2012 onward",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  modern_icd10 <- primary[
    index_date >= as.Date("2015-10-01")
  ]

  res <- run_ato_generic(
    modern_icd10,
    "Modern-era primary cohort: ICD-10 era index",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  # ----------------------------------------------------------
  # 4. NEGATIVE-CONTROL OUTCOMES
  # ----------------------------------------------------------
  run_nco_if_feasible <- function(
      tbl,
      label,
      min_total_events = 20L
  ) {

    d <- read_table_cols(
      tbl,
      c(
        "nco_event",
        "nco_followup_days"
      )
    )

    d[, trt := as.integer(exposure == "5ARI")]

    n_events <- sum(d$nco_event, na.rm = TRUE)
    events_treated <- sum(
      d$nco_event[d$trt == 1],
      na.rm = TRUE
    )
    events_control <- sum(
      d$nco_event[d$trt == 0],
      na.rm = TRUE
    )

    if (
      n_events < min_total_events ||
      events_treated < 3L ||
      events_control < 3L
    ) {
      warning(
        label,
        ": too sparse for a stable negative-control HR; ",
        "counts will remain available from Stage 14."
      )
      return(NULL)
    }

    run_ato_generic(
      d,
      label,
      "nco_event",
      "nco_followup_days"
    )
  }

  res <- run_nco_if_feasible(
    "analysis_cohort_nco_cranial_nerve",
    "Negative control: benign cranial nerve neoplasm (D33.3 / ICD-9 225.1)"
  )
  if (!is.null(res)) store_result(res)

  res <- run_nco_if_feasible(
    "analysis_cohort_nco_spinal_cord",
    "Negative control: benign spinal cord neoplasm (D33.4 / ICD-9 225.3)"
  )
  if (!is.null(res)) store_result(res)

  # ----------------------------------------------------------
  # 5. EHR-CONTINUITY IPCW SENSITIVITY
  # ----------------------------------------------------------
  #
  # We explicitly model loss of observable EHR follow-up in 180-d
  # intervals. The denominator model includes:
  #   - treatment
  #   - baseline confounders
  #   - baseline utilization
  #   - time since risk-set entry
  #   - healthcare utilization during the immediately PRIOR 180-d interval
  #
  # The numerator model includes treatment + time only.
  #
  # Weights at each interval are based on the probability of remaining
  # uncensored THROUGH PRIOR intervals, so current-interval censoring is
  # not used to weight its own at-risk contribution.
  #
  fit_large_binomial <- function(
      formula,
      data,
      label
  ) {

    if (requireNamespace("speedglm", quietly = TRUE)) {

      fit <- try(
        speedglm::speedglm(
          formula,
          data = data,
          family = binomial()
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        return(fit)
      }

      warning(
        label,
        ": speedglm failed; falling back to glm()."
      )
    }

    glm(
      formula,
      data = data,
      family = binomial(),
      control = glm.control(maxit = 50)
    )
  }

  predict_prob <- function(fit, newdata) {

    p <- as.numeric(
      predict(
        fit,
        newdata = newdata,
        type = "response"
      )
    )

    pmin(
      pmax(p, 1e-6),
      1 - 1e-6
    )
  }

  extract_ipcw_hr <- function(
      pp,
      weight_col,
      label
  ) {

    fit <- coxph(
      Surv(
        tstart_y,
        tstop_y,
        outcome_event
      ) ~ trt + cluster(patient_id),
      data = pp,
      weights = pp[[weight_col]],
      ties = "efron"
    )

    h <- extract_hr(fit)

    data.table(
      analysis = label,
      n_patients = uniqueN(pp$patient_id),
      person_period_rows = nrow(pp),
      events = sum(pp$outcome_event),
      HR = h$HR,
      CI_low = h$CI_low,
      CI_high = h$CI_high,
      p = h$p,
      weight = weight_col
    )
  }

  cat("\n============================================================\n")
  cat("Building 180-d person-period data for IPCW sensitivity...\n")
  cat("============================================================\n")

  # Re-estimate the finalized primary overlap weights on the exact IPCW
  # analysis population.
  primary_ipcw <- copy(primary)

  primary_ipcw[
    ,
    `:=`(
      last_encounter_date = as.Date(last_encounter_date),
      landmark_date = as.Date(landmark_date),
      incident_meningioma_date_strict =
        as.Date(incident_meningioma_date_strict)
    )
  ]

  primary_ipcw <- primary_ipcw[
    is.finite(followup_days_strict) &
    followup_days_strict > 0
  ]

  wt <- estimate_overlap_weights(
    primary_ipcw,
    label = "Primary cohort for EHR-continuity IPCW"
  )

  b <- wt$data

  data_cutoff <- as.Date(
    dbGetQuery(
      con,
      "SELECT MAX(last_encounter_date) AS x
       FROM candidate_observation;"
    )$x[1]
  )

  admin_threshold <- data_cutoff - ADMIN_GRACE_DAYS

  b[
    ,
    n_intervals := pmax(
      1L,
      as.integer(
        ceiling(
          as.numeric(followup_days_strict) /
          IPCW_INTERVAL_DAYS
        )
      )
    )
  ]

  row_index <- rep(
    seq_len(nrow(b)),
    times = b$n_intervals
  )

  interval_index <- sequence(
    b$n_intervals
  ) - 1L

  pp <- b[row_index]
  pp[, interval_index := interval_index]

  pp[
    ,
    `:=`(
      tstart_days =
        interval_index * IPCW_INTERVAL_DAYS,

      tstop_days = pmin(
        (interval_index + 1L) * IPCW_INTERVAL_DAYS,
        as.numeric(followup_days_strict)
      )
    )
  ]

  pp[
    ,
    outcome_event := as.integer(
      meningioma_event_strict == 1L &
      as.numeric(followup_days_strict) > tstart_days &
      as.numeric(followup_days_strict) <= tstop_days
    )
  ]

  pp[
    ,
    final_interval := as.integer(
      tstop_days >= as.numeric(followup_days_strict)
    )
  ]

  pp[
    ,
    censor_event := as.integer(
      final_interval == 1L &
      outcome_event == 0L &
      last_encounter_date < admin_threshold
    )
  ]

  pp[
    ,
    administrative_end := as.integer(
      final_interval == 1L &
      outcome_event == 0L &
      censor_event == 0L
    )
  ]

  pp[
    ,
    prev_bin := interval_index - 1L
  ]

  util <- as.data.table(
    dbGetQuery(
      con,
"
SELECT
  patient_id,
  interval_bin,
  encounters_total,
  encounter_days,
  ambulatory_encounters,
  ed_encounters,
  inpatient_encounters
FROM method_primary_180d_util;
"
    )
  )

  setnames(
    util,
    c(
      "encounters_total",
      "encounter_days",
      "ambulatory_encounters",
      "ed_encounters",
      "inpatient_encounters"
    ),
    c(
      "prev_encounters_total",
      "prev_encounter_days",
      "prev_ambulatory",
      "prev_ed",
      "prev_inpatient"
    )
  )

  pp <- merge(
    pp,
    util,
    by.x = c("patient_id", "prev_bin"),
    by.y = c("patient_id", "interval_bin"),
    all.x = TRUE,
    sort = FALSE
  )

  # PATCH: enforce one row per patient / interval and chronological order.
  # This protects the subsequent within-patient cumulative IPCW calculation
  # from implicit ordering assumptions after the utilization merge.
  if (anyDuplicated(pp[, .(patient_id, interval_index)])) {
    stop("IPCW person-period merge produced duplicate patient/interval rows.")
  }
  data.table::setorder(pp, patient_id, interval_index)

  util_vars <- c(
    "prev_encounters_total",
    "prev_encounter_days",
    "prev_ambulatory",
    "prev_ed",
    "prev_inpatient"
  )

  for (v in util_vars) {
    set(
      pp,
      i = which(is.na(pp[[v]])),
      j = v,
      value = 0
    )

    pp[
      ,
      paste0("log_", v) :=
        log1p(as.numeric(get(v)))
    ]
  }

  pp[
    ,
    interval_t :=
      interval_index /
      max(1, max(interval_index, na.rm = TRUE))
  ]

  # Dynamically keep only varying baseline terms.
  factor_terms <- FACTOR_CANDIDATES[
    vapply(
      FACTOR_CANDIDATES,
      function(v) {
        v %in% names(pp) &&
        n_observed_unique(pp[[v]]) >= 2L
      },
      logical(1)
    )
  ]

  binary_terms <- BINARY_CANDIDATES[
    vapply(
      BINARY_CANDIDATES,
      function(v) {
        v %in% names(pp) &&
        n_observed_unique(pp[[v]]) >= 2L
      },
      logical(1)
    )
  ]

  baseline_util_terms <- UTILIZATION_CANDIDATES[
    vapply(
      UTILIZATION_CANDIDATES,
      function(v) {
        v %in% names(pp) &&
        n_observed_unique(pp[[v]]) >= 2L
      },
      logical(1)
    )
  ]

  prior_util_terms <- paste0(
    "log_",
    util_vars
  )

  den_terms <- c(
    "trt",
    "interval_t",
    "I(interval_t^2)",
    "I(interval_t^3)",
    "splines::ns(age_at_index, df=5)",
    "splines::ns(index_year, df=5)",
    factor_terms,
    binary_terms,
    baseline_util_terms,
    prior_util_terms
  )

  den_formula <- as.formula(
    paste(
      "censor_event ~",
      paste(unique(den_terms), collapse = " + ")
    )
  )

  num_formula <- censor_event ~
    trt +
    interval_t +
    I(interval_t^2) +
    I(interval_t^3)

  censor_events <- sum(pp$censor_event)
  if (censor_events < 100L) {
    stop(
      "Too few loss-to-observation censoring events for IPCW: ",
      censor_events
    )
  }

  cat(
    "Person-period rows: ",
    format(nrow(pp), big.mark = ","),
    "\n",
    sep = ""
  )

  cat(
    "Loss-to-observation censoring events: ",
    format(censor_events, big.mark = ","),
    "\n",
    sep = ""
  )

  cat("\nFitting denominator censoring model...\n")
  den_fit <- fit_large_binomial(
    den_formula,
    pp,
    "IPCW denominator"
  )

  cat("Fitting numerator censoring model...\n")
  num_fit <- fit_large_binomial(
    num_formula,
    pp,
    "IPCW numerator"
  )

  pp[, p_censor_den := predict_prob(den_fit, pp)]
  pp[, p_censor_num := predict_prob(num_fit, pp)]

  pp[
    ,
    log_keep_ratio :=
      log1p(-p_censor_num) -
      log1p(-p_censor_den)
  ]

  data.table::setorder(pp, patient_id, interval_index)

  # Weight each interval by probability of having remained uncensored
  # through all PRIOR completed intervals.
  pp[
    ,
    log_ipcw_raw := shift(
      cumsum(log_keep_ratio),
      n = 1L,
      fill = 0
    ),
    by = patient_id
  ]

  pp[
    ,
    ipcw_raw := exp(
      pmin(
        pmax(log_ipcw_raw, -20),
        20
      )
    )
  ]

  q_1_99 <- quantile(
    pp$ipcw_raw,
    probs = c(0.01, 0.99),
    na.rm = TRUE
  )

  q_025_975 <- quantile(
    pp$ipcw_raw,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )

  pp[
    ,
    ipcw_1_99 := pmin(
      pmax(ipcw_raw, q_1_99[1]),
      q_1_99[2]
    )
  ]

  pp[
    ,
    ipcw_025_975 := pmin(
      pmax(ipcw_raw, q_025_975[1]),
      q_025_975[2]
    )
  ]

  pp[
    ,
    `:=`(
      final_w_raw = ow * ipcw_raw,
      final_w_1_99 = ow * ipcw_1_99,
      final_w_025_975 = ow * ipcw_025_975,
      tstart_y = tstart_days / 365.25,
      tstop_y = tstop_days / 365.25
    )
  ]

  ipcw_results <- rbindlist(
    list(
      extract_ipcw_hr(
        pp,
        "final_w_raw",
        "Primary + EHR-continuity IPCW: untruncated"
      ),
      extract_ipcw_hr(
        pp,
        "final_w_1_99",
        "Primary + EHR-continuity IPCW: censoring weights truncated 1st/99th pct"
      ),
      extract_ipcw_hr(
        pp,
        "final_w_025_975",
        "Primary + EHR-continuity IPCW: censoring weights truncated 2.5th/97.5th pct"
      )
    ),
    fill = TRUE
  )

  ipcw_diag <- data.table(
    metric = c(
      "patients",
      "person_period_rows",
      "outcome_events",
      "loss_to_observation_censor_events",
      "administrative_end_rows",
      "ipcw_raw_p01",
      "ipcw_raw_p50",
      "ipcw_raw_p99",
      "ipcw_raw_max",
      "combined_weight_1_99_p01",
      "combined_weight_1_99_p50",
      "combined_weight_1_99_p99",
      "combined_weight_1_99_max",
      "predicted_censor_den_p01",
      "predicted_censor_den_p50",
      "predicted_censor_den_p99"
    ),
    value = c(
      uniqueN(pp$patient_id),
      nrow(pp),
      sum(pp$outcome_event),
      sum(pp$censor_event),
      sum(pp$administrative_end),
      quantile(pp$ipcw_raw, 0.01, na.rm = TRUE),
      quantile(pp$ipcw_raw, 0.50, na.rm = TRUE),
      quantile(pp$ipcw_raw, 0.99, na.rm = TRUE),
      max(pp$ipcw_raw, na.rm = TRUE),
      quantile(pp$final_w_1_99, 0.01, na.rm = TRUE),
      quantile(pp$final_w_1_99, 0.50, na.rm = TRUE),
      quantile(pp$final_w_1_99, 0.99, na.rm = TRUE),
      max(pp$final_w_1_99, na.rm = TRUE),
      quantile(pp$p_censor_den, 0.01, na.rm = TRUE),
      quantile(pp$p_censor_den, 0.50, na.rm = TRUE),
      quantile(pp$p_censor_den, 0.99, na.rm = TRUE)
    )
  )

  # ----------------------------------------------------------
  # 6. EXPLORATORY-FAMILY FDR
  # ----------------------------------------------------------
  stage12_file <- file.path(
    SIGNAL_DIR,
    "12_signal_model_results.csv"
  )

  if (file.exists(stage12_file)) {

    exploratory <- fread(stage12_file)

    if (!"weighted_p" %in% names(exploratory)) {
      stop(
        "12_signal_model_results.csv does not contain weighted_p."
      )
    }

    exploratory[
      ,
      q_BH := p.adjust(
        weighted_p,
        method = "BH"
      )
    ]

    exploratory[
      ,
      `:=`(
        raw_p_lt_0_05 =
          as.integer(weighted_p < 0.05),
        q_BH_lt_0_05 =
          as.integer(q_BH < 0.05)
      )
    ]

    fwrite(
      exploratory,
      file.path(
        OUT_DIR,
        "15_stage12_exploratory_FDR.csv"
      )
    )
  } else {

    warning(
      "Stage-12 model file not found; FDR table not produced."
    )
  }

  # ----------------------------------------------------------
  # 7. SAVE STANDARD ATO RESULTS
  # ----------------------------------------------------------
  models <- rbindlist(
    model_rows,
    fill = TRUE
  )

  balances <- rbindlist(
    balance_rows,
    fill = TRUE
  )

  diagnostics <- rbindlist(
    diag_rows,
    fill = TRUE
  )

  incidence <- rbindlist(
    incidence_rows,
    fill = TRUE
  )

  fwrite(
    models,
    file.path(
      OUT_DIR,
      "15_methodologic_model_results.csv"
    )
  )

  fwrite(
    balances,
    file.path(
      OUT_DIR,
      "15_methodologic_balance_all.csv"
    )
  )

  fwrite(
    balance_summary(balances),
    file.path(
      OUT_DIR,
      "15_methodologic_balance_summary.csv"
    )
  )

  fwrite(
    diagnostics,
    file.path(
      OUT_DIR,
      "15_methodologic_predictor_diagnostics.csv"
    )
  )

  fwrite(
    incidence,
    file.path(
      OUT_DIR,
      "15_methodologic_incidence_rates.csv"
    )
  )

  if (length(ph_rows) > 0L) {
    fwrite(
      rbindlist(ph_rows, fill = TRUE),
      file.path(
        OUT_DIR,
        "15_methodologic_ph_diagnostics.csv"
      )
    )
  }

  fwrite(
    ipcw_results,
    file.path(
      OUT_DIR,
      "15_ipcw_results.csv"
    )
  )

  fwrite(
    ipcw_diag,
    file.path(
      OUT_DIR,
      "15_ipcw_diagnostics.csv"
    )
  )

  # ----------------------------------------------------------
  # 8. COMPACT REVIEW FILE
  # ----------------------------------------------------------
  review <- rbindlist(
    list(
      models[
        ,
        .(
          analysis,
          n,
          events,
          weighted_HR,
          weighted_CI_low,
          weighted_CI_high,
          weighted_p
        )
      ],

      ipcw_results[
        ,
        .(
          analysis,
          n = n_patients,
          events,
          weighted_HR = HR,
          weighted_CI_low = CI_low,
          weighted_CI_high = CI_high,
          weighted_p = p
        )
      ]
    ),
    fill = TRUE
  )

  fwrite(
    review,
    file.path(
      OUT_DIR,
      "15_REVIEW_FIRST.csv"
    )
  )

  cat("\n================ STANDARD TARGETED MODELS ============\n")
  print(models)

  cat("\n================ IPCW SENSITIVITY ====================\n")
  print(ipcw_results)

  cat("\n================ BALANCE SUMMARY =====================\n")
  print(balance_summary(balances))

  cat("\nStage 15 complete.\n")
  cat("Analysis directory:", OUT_DIR, "\n")
  cat(
    "\nFIRST FILE TO SEND BACK:\n",
    file.path(OUT_DIR, "15_REVIEW_FIRST.csv"),
    "\n",
    sep = ""
  )
}

run_stage_15()
