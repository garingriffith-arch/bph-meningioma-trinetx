# NOTE: This stage is wrapped in a function so on.exit() is tied to
# the stage function rather than source()/eval(). This prevents the
# DuckDB connection from being closed immediately when sourced.

run_stage_10 <- function() {
  # ============================================================
  # 10_run_final_sensitivity_models.R
  #
  # Runs the final literature-informed analysis package:
  #   - sustained 1-y treatment
  #   - strict sustained 1-y treatment
  #   - persistent 3-y treatment
  #   - high-specificity meningioma outcome
  #   - maximally permissive 1-y-lag model
  #   - permissive-BPH 1-y-lag model
  #   - maximally permissive zero-lag model
  #   - finasteride vs dutasteride
  #   - finasteride 5 mg vs dutasteride 0.5 mg
  #   - interval-specific latency analyses
  #   - registry-confirmed outcomes when available
  #   - treatment-requiring outcomes when available
  #   - death competing-risk Fine-Gray sensitivity
  #   - measured-BMI sensitivity when sufficiently complete
  #   - days-supply sustained sensitivity when feasible
  #
  # Primary estimand for binary treatment contrasts:
  #   ATO via logistic propensity-score overlap weighting.
  # ============================================================

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
  OUT_DIR <- file.path(WORK_DIR, "analysis_final")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  source(file.path(R_DIR, "final_analysis_helpers.R"))

  N_THREADS <- min(24L, max(4L, parallel::detectCores(logical = TRUE) - 6L))
  MEMORY_LIMIT <- "48GB"

  con <- dbConnect(duckdb::duckdb(), dbdir = DB_PATH, read_only = FALSE)
  on.exit({
    try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }, add = TRUE)

  dbExecute(con, sprintf("PRAGMA threads=%d;", N_THREADS))
  dbExecute(con, sprintf("PRAGMA memory_limit='%s';", MEMORY_LIMIT))

  tbls <- dbListTables(con)

  required <- c(
    "analysis_cohort_primary_revised",
    "analysis_cohort_sustained_1y",
    "analysis_cohort_sustained_1y_strict",
    "analysis_cohort_persistent_3y",
    "analysis_cohort_specific_meningioma",
    "analysis_cohort_competing_risk",
    "analysis_cohort_permissive_max",
    "analysis_cohort_permissive_bph",
    "analysis_cohort_permissive_zero_lag",
    "analysis_cohort_latency_clean"
  )

  missing <- setdiff(required, tbls)
  if (length(missing) > 0L) {
    stop(
      "Run Stages 07 and 08 first. Missing: ",
      paste(missing, collapse=", ")
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

  read_table_cols <- function(tbl, extra_cols) {
    fields <- dbListFields(con, tbl)
    use <- unique(c(common_cols, extra_cols))
    use <- use[use %in% fields]

    as.data.table(
      dbGetQuery(
        con,
        paste0(
          "SELECT ",
          paste(use, collapse=", "),
          " FROM ",
          tbl,
          ";"
        )
      )
    )
  }

  analysis_objects <- list()
  model_rows <- list()
  balance_rows <- list()
  diag_rows <- list()
  incidence_rows <- list()
  ph_rows <- list()

  store_result <- function(x) {
    analysis_objects[[length(analysis_objects) + 1L]] <<- x
    model_rows[[length(model_rows) + 1L]] <<- x$result
    balance_rows[[length(balance_rows) + 1L]] <<- x$balance
    diag_rows[[length(diag_rows) + 1L]] <<- x$diagnostics
    incidence_rows[[length(incidence_rows) + 1L]] <<- x$incidence
    if (!is.null(x$zph)) {
      ph_rows[[length(ph_rows) + 1L]] <<- x$zph
    }
  }

  # ------------------------------------------------------------
  # 1. SUSTAINED 1-Y TREATMENT
  # ------------------------------------------------------------
  s1 <- read_table_cols(
    "analysis_cohort_sustained_1y",
    c("meningioma_event_strict", "followup_days_strict")
  )
  s1[, trt := as.integer(exposure=="5ARI")]

  res <- run_ato_generic(
    s1,
    "Sustained treatment: repeat same-class day 90-365, no crossover",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  # ------------------------------------------------------------
  # 2. STRICTER 1-Y SUSTAINED
  # ------------------------------------------------------------
  s1s <- read_table_cols(
    "analysis_cohort_sustained_1y_strict",
    c("meningioma_event_strict", "followup_days_strict")
  )
  s1s[, trt := as.integer(exposure=="5ARI")]

  res <- run_ato_generic(
    s1s,
    "Sustained treatment strict: repeat evidence in both 90-180 and 181-365 d",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  # ------------------------------------------------------------
  # 3. 3-Y PERSISTENT TREATMENT
  # ------------------------------------------------------------
  p3 <- read_table_cols(
    "analysis_cohort_persistent_3y",
    c("meningioma_event_strict", "followup_days_strict")
  )
  p3[, trt := as.integer(exposure=="5ARI")]

  if (
    nrow(p3[trt==1]) >= 500 &&
    nrow(p3[trt==0]) >= 500 &&
    sum(p3$meningioma_event_strict, na.rm=TRUE) >= 40
  ) {
    res <- run_ato_generic(
      p3,
      "Persistent treatment: same-class evidence in years 1, 2, and 3; 3-y landmark",
      "meningioma_event_strict",
      "followup_days_strict"
    )
    store_result(res)
  } else {
    warning("3-y persistence cohort too sparse for stable modeling; skipped.")
  }

  # ------------------------------------------------------------
  # 4. HIGH-SPECIFICITY MENINGIOMA
  # ------------------------------------------------------------
  sp <- read_table_cols(
    "analysis_cohort_specific_meningioma",
    c(
      "meningioma_event_specific",
      "followup_days_specific",
      "meningioma_event_d320",
      "followup_days_d320"
    )
  )
  sp[, trt := as.integer(exposure=="5ARI")]

  res <- run_ato_generic(
    sp,
    "High-specificity intracranial meningioma: D32.0 + ICD-9 225.2",
    "meningioma_event_specific",
    "followup_days_specific"
  )
  store_result(res)

  # Pure D32.0 analysis is restricted to ICD-10-era initiators to avoid giving
  # the 5ARI/control groups unequal historical opportunity for D32.0 coding.
  d320 <- sp[index_date >= as.Date("2015-10-01")]
  if (sum(d320$meningioma_event_d320, na.rm=TRUE) >= 40) {
    res <- run_ato_generic(
      d320,
      "ICD-10-era D32.0-only meningioma",
      "meningioma_event_d320",
      "followup_days_d320"
    )
    store_result(res)
  }

  # ------------------------------------------------------------
  # 5. MAXIMALLY PERMISSIVE 1-Y-LAG MODEL
  # ------------------------------------------------------------
  pm <- read_table_cols(
    "analysis_cohort_permissive_max",
    c(
      "prior_prostate_cancer",
      "bph_near_index",
      "known_finasteride_1mg",
      "meningioma_event_strict",
      "followup_days_strict"
    )
  )
  pm[, trt := as.integer(exposure=="5ARI")]

  res <- run_ato_generic(
    pm,
    "Maximally permissive: any 5ARI-containing first regimen vs uro-alpha, 1-y lag",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  # ------------------------------------------------------------
  # 6. PERMISSIVE BPH 1-Y-LAG MODEL
  # ------------------------------------------------------------
  pb <- read_table_cols(
    "analysis_cohort_permissive_bph",
    c(
      "known_finasteride_1mg",
      "meningioma_event_strict",
      "followup_days_strict"
    )
  )
  pb[, trt := as.integer(exposure=="5ARI")]

  res <- run_ato_generic(
    pb,
    "Permissive BPH: combination starts/crossover/1-mg allowed, 1-y lag",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  # ------------------------------------------------------------
  # 7. MAXIMALLY PERMISSIVE ZERO-LAG MODEL
  # ------------------------------------------------------------
  pz <- read_table_cols(
    "analysis_cohort_permissive_zero_lag",
    c(
      "prior_prostate_cancer",
      "bph_near_index",
      "known_finasteride_1mg",
      "meningioma_event_strict",
      "followup_days_strict"
    )
  )
  pz[, trt := as.integer(exposure=="5ARI")]

  res <- run_ato_generic(
    pz,
    "Maximally permissive zero-lag",
    "meningioma_event_strict",
    "followup_days_strict"
  )
  store_result(res)

  # ------------------------------------------------------------
  # 8. FINASTERIDE VS DUTASTERIDE HEAD-TO-HEAD
  # ------------------------------------------------------------
  primary <- read_table_cols(
    "analysis_cohort_primary_revised",
    c("meningioma_event_strict", "followup_days_strict")
  )

  fd <- primary[
    exposure=="5ARI" &
    index_drug %in% c("finasteride","dutasteride")
  ]
  fd[, trt := as.integer(index_drug=="finasteride")]

  if (
    nrow(fd[trt==1]) >= 500 &&
    nrow(fd[trt==0]) >= 500 &&
    sum(fd$meningioma_event_strict, na.rm=TRUE) >= 40
  ) {
    res <- run_ato_generic(
      fd,
      "Head-to-head: finasteride vs dutasteride",
      "meningioma_event_strict",
      "followup_days_strict"
    )
    store_result(res)
  }

  fd_dose <- primary[
    exposure=="5ARI" &
    (
      (
        index_drug=="finasteride" &
        !is.na(index_finasteride_strength) &
        toupper(trimws(index_finasteride_strength))=="5 MG"
      )
      |
      (
        index_drug=="dutasteride" &
        !is.na(index_dutasteride_strength) &
        toupper(trimws(index_dutasteride_strength))=="0.5 MG"
      )
    )
  ]
  fd_dose[, trt := as.integer(index_drug=="finasteride")]

  if (
    nrow(fd_dose[trt==1]) >= 500 &&
    nrow(fd_dose[trt==0]) >= 500 &&
    sum(fd_dose$meningioma_event_strict, na.rm=TRUE) >= 30
  ) {
    res <- run_ato_generic(
      fd_dose,
      "Head-to-head exact BPH doses: finasteride 5 mg vs dutasteride 0.5 mg",
      "meningioma_event_strict",
      "followup_days_strict"
    )
    store_result(res)
  }

  # ------------------------------------------------------------
  # 9. INDEPENDENT REGISTRY ENDPOINTS
  # ------------------------------------------------------------
  if ("analysis_cohort_registry_meningioma" %in% tbls) {

    rg <- read_table_cols(
      "analysis_cohort_registry_meningioma",
      c(
        "registry_meningioma_event",
        "registry_followup_days",
        "cancer_registry_meningioma_event",
        "cancer_registry_followup_days"
      )
    )
    rg[, trt := as.integer(exposure=="5ARI")]

    if (sum(rg$registry_meningioma_event, na.rm=TRUE) >= 30) {
      res <- run_ato_generic(
        rg,
        "ICD-O tumor-table meningioma",
        "registry_meningioma_event",
        "registry_followup_days"
      )
      store_result(res)
    }

    if (sum(rg$cancer_registry_meningioma_event, na.rm=TRUE) >= 30) {
      res <- run_ato_generic(
        rg,
        "Cancer-registry-source ICD-O meningioma",
        "cancer_registry_meningioma_event",
        "cancer_registry_followup_days"
      )
      store_result(res)
    }
  }

  # ------------------------------------------------------------
  # 10. CLINICALLY TREATED MENINGIOMA
  # ------------------------------------------------------------
  if ("analysis_cohort_treated_meningioma" %in% tbls) {

    tx <- read_table_cols(
      "analysis_cohort_treated_meningioma",
      c("treated_meningioma_event", "treated_meningioma_followup_days")
    )
    tx[, trt := as.integer(exposure=="5ARI")]

    if (sum(tx$treated_meningioma_event, na.rm=TRUE) >= 30) {
      res <- run_ato_generic(
        tx,
        "Clinically treated meningioma: resection or cranial SRS within 180 d",
        "treated_meningioma_event",
        "treated_meningioma_followup_days"
      )
      store_result(res)
    }
  }

  if ("analysis_cohort_oncology_treated_registry" %in% tbls) {

    ot <- read_table_cols(
      "analysis_cohort_oncology_treated_registry",
      c(
        "oncology_treated_registry_event",
        "oncology_treated_registry_followup_days"
      )
    )
    ot[, trt := as.integer(exposure=="5ARI")]

    if (sum(ot$oncology_treated_registry_event, na.rm=TRUE) >= 30) {
      res <- run_ato_generic(
        ot,
        "Oncology-treatment-table ICD-O meningioma",
        "oncology_treated_registry_event",
        "oncology_treated_registry_followup_days"
      )
      store_result(res)
    }
  }

  # ------------------------------------------------------------
  # 11. MEASURED BMI SENSITIVITY
  # ------------------------------------------------------------
  if ("analysis_cohort_primary_bmi" %in% tbls) {

    bmi <- read_table_cols(
      "analysis_cohort_primary_bmi",
      c(
        "bmi_preindex",
        "bmi_missing",
        "meningioma_event_strict",
        "followup_days_strict"
      )
    )
    bmi[, trt := as.integer(exposure=="5ARI")]

    bmi_qc <- bmi[, .(
      observed_pct=100*mean(!is.na(bmi_preindex))
    ), by=trt]

    fwrite(
      bmi_qc,
      file.path(OUT_DIR, "10_bmi_model_feasibility.csv")
    )

    # Exploratory missing-indicator sensitivity. We only run it when at least
    # half of each arm has measured BMI; otherwise coded obesity remains the
    # more defensible primary covariate.
    if (all(bmi_qc$observed_pct >= 50)) {
      med_bmi <- median(bmi$bmi_preindex, na.rm=TRUE)
      bmi[, bmi_missing := as.integer(is.na(bmi_preindex))]
      bmi[, bmi_imputed := fifelse(
        is.na(bmi_preindex),
        med_bmi,
        as.numeric(bmi_preindex)
      )]

      res <- run_ato_generic(
        bmi,
        "Measured-BMI sensitivity (median imputation + missingness indicator)",
        "meningioma_event_strict",
        "followup_days_strict",
        extra_numeric="bmi_imputed"
      )
      store_result(res)
    }
  }

  # ------------------------------------------------------------
  # 12. DAYS-SUPPLY EXPLORATORY SUSTAINED COHORT
  # ------------------------------------------------------------
  if ("analysis_cohort_days_supply_exploratory" %in% tbls) {

    ds <- read_table_cols(
      "analysis_cohort_days_supply_exploratory",
      c(
        "meningioma_event_strict",
        "followup_days_strict",
        "approx_days_supply_sum_365"
      )
    )
    ds[, trt := as.integer(exposure=="5ARI")]

    if (
      nrow(ds[trt==1]) >= 1000 &&
      nrow(ds[trt==0]) >= 1000 &&
      sum(ds$meningioma_event_strict, na.rm=TRUE) >= 40
    ) {
      res <- run_ato_generic(
        ds,
        "Exploratory days-supply sustained cohort: >=180 approximate days in first year",
        "meningioma_event_strict",
        "followup_days_strict"
      )
      store_result(res)
    }
  }

  # ------------------------------------------------------------
  # 13. LATENCY / TIME-SINCE-INITIATION ANALYSIS
  # ------------------------------------------------------------
  lat <- read_table_cols(
    "analysis_cohort_latency_clean",
    c("meningioma_event_strict", "followup_days_strict")
  )
  lat[, trt := as.integer(exposure=="5ARI")]
  lat[, total_time_y := as.numeric(followup_days_strict)/365.25]
  lat[, event := as.integer(meningioma_event_strict)]

  wt_lat <- estimate_overlap_weights(
    lat,
    label="Latency clean cohort baseline weights"
  )
  latw <- wt_lat$data

  intervals <- data.table(
    interval=c("0-<1 y","1-<2 y","2-<5 y",">=5 y"),
    start=c(0,1,2,5),
    end=c(1,2,5,Inf)
  )

  latency_rows <- list()

  for (i in seq_len(nrow(intervals))) {

    a <- intervals$start[i]
    b <- intervals$end[i]

    x <- latw[total_time_y > a]
    if (nrow(x) == 0L) next

    if (is.finite(b)) {
      x[, interval_time := pmin(total_time_y, b) - a]
      x[, interval_event := as.integer(
        event==1 &
        total_time_y > a &
        total_time_y <= b
      )]
    } else {
      x[, interval_time := total_time_y - a]
      x[, interval_event := as.integer(
        event==1 &
        total_time_y > a
      )]
    }

    x <- x[interval_time > 0]

    if (
      sum(x$interval_event) < 10 ||
      length(unique(x$trt)) < 2L
    ) {
      latency_rows[[length(latency_rows)+1L]] <- data.table(
        interval=intervals$interval[i],
        n=nrow(x),
        events=sum(x$interval_event),
        HR=NA_real_,
        CI_low=NA_real_,
        CI_high=NA_real_,
        p=NA_real_
      )
      next
    }

    fit <- coxph(
      Surv(interval_time, interval_event) ~ trt,
      data=x,
      weights=ow,
      robust=TRUE,
      ties="efron"
    )

    h <- extract_hr(fit)

    latency_rows[[length(latency_rows)+1L]] <- data.table(
      interval=intervals$interval[i],
      n=nrow(x),
      n_5ari=sum(x$trt==1),
      n_uro_alpha=sum(x$trt==0),
      events=sum(x$interval_event),
      events_5ari=sum(x$interval_event[x$trt==1]),
      events_uro_alpha=sum(x$interval_event[x$trt==0]),
      HR=h$HR,
      CI_low=h$CI_low,
      CI_high=h$CI_high,
      p=h$p
    )
  }

  latency_results <- rbindlist(latency_rows, fill=TRUE)
  fwrite(
    latency_results,
    file.path(OUT_DIR, "10_latency_interval_results.csv")
  )

  fwrite(
    wt_lat$balance,
    file.path(OUT_DIR, "10_latency_balance.csv")
  )

  # ------------------------------------------------------------
  # 14. COMPETING-RISK FINE-GRAY
  # ------------------------------------------------------------
  comp <- read_table_cols(
    "analysis_cohort_competing_risk",
    c("competing_status", "competing_time_days")
  )
  comp[, trt := as.integer(exposure=="5ARI")]
  comp[, competing_time_y := as.numeric(competing_time_days)/365.25]
  comp <- comp[
    is.finite(competing_time_y) &
    competing_time_y > 0 &
    competing_status %in% c(0,1,2)
  ]

  competing_result <- NULL

  if (
    sum(comp$competing_status==1) >= 30 &&
    sum(comp$competing_status==2) >= 30
  ) {
    wt_comp <- estimate_overlap_weights(
      comp,
      label="Competing-risk baseline weights"
    )
    cw <- wt_comp$data

    cw[, status_factor := factor(
      competing_status,
      levels=c(0,1,2),
      labels=c("censor","meningioma","death")
    )]
    cw[, fg_id := .I]

    fg_try <- try(
      finegray(
        Surv(competing_time_y, status_factor) ~ trt,
        data=cw,
        etype="meningioma",
        id=fg_id,
        weights=ow
      ),
      silent=TRUE
    )

    if (!inherits(fg_try, "try-error")) {

      if ("fg_id" %in% names(fg_try)) {
        fg_fit <- coxph(
          Surv(fgstart, fgstop, fgstatus) ~ trt + cluster(fg_id),
          data=fg_try,
          weights=fgwt,
          ties="efron"
        )
      } else {
        warning(
          "finegray() did not retain fg_id; fitting robust weighted Fine-Gray ",
          "without explicit subject clustering."
        )
        fg_fit <- coxph(
          Surv(fgstart, fgstop, fgstatus) ~ trt,
          data=fg_try,
          weights=fgwt,
          robust=TRUE,
          ties="efron"
        )
      }

      h <- extract_hr(fg_fit)

      competing_result <- data.table(
        analysis="Fine-Gray competing-risk sensitivity: death competing event",
        n=nrow(cw),
        meningioma_events=sum(cw$competing_status==1),
        deaths=sum(cw$competing_status==2),
        subdistribution_HR=h$HR,
        CI_low=h$CI_low,
        CI_high=h$CI_high,
        p=h$p
      )

      fwrite(
        competing_result,
        file.path(OUT_DIR, "10_competing_risk_result.csv")
      )
    } else {
      writeLines(
        as.character(fg_try),
        file.path(OUT_DIR, "10_competing_risk_error.txt")
      )
    }
  }

  # ------------------------------------------------------------
  # 15. SAVE MASTER OUTPUTS
  # ------------------------------------------------------------
  models <- rbindlist(model_rows, fill=TRUE)
  balances <- rbindlist(balance_rows, fill=TRUE)
  diagnostics <- rbindlist(diag_rows, fill=TRUE)
  incidence <- rbindlist(incidence_rows, fill=TRUE)
  ph <- if (length(ph_rows) > 0L) {
    rbindlist(ph_rows, fill=TRUE)
  } else {
    data.table()
  }

  fwrite(models, file.path(OUT_DIR, "10_final_model_results.csv"))
  fwrite(balances, file.path(OUT_DIR, "10_final_balance_all.csv"))
  fwrite(
    balance_summary(balances),
    file.path(OUT_DIR, "10_final_balance_summary.csv")
  )
  fwrite(
    diagnostics,
    file.path(OUT_DIR, "10_final_predictor_diagnostics.csv")
  )
  fwrite(
    incidence,
    file.path(OUT_DIR, "10_final_incidence_rates.csv")
  )
  if (nrow(ph) > 0L) {
    fwrite(ph, file.path(OUT_DIR, "10_final_ph_diagnostics.csv"))
  }

  cat("\n================ FINAL MODEL RESULTS ================\n")
  print(models)

  cat("\n================ LATENCY RESULTS ====================\n")
  print(latency_results)

  if (!is.null(competing_result)) {
    cat("\n================ COMPETING RISK =====================\n")
    print(competing_result)
  }

  cat("\nFinal sensitivity analysis complete.\n")
  cat("Outputs:", OUT_DIR, "\n")

  dbDisconnect(con, shutdown = TRUE)
}

run_stage_10()
