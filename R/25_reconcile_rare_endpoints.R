# Stage 25: Final-cohort rare-endpoint reconciliation.
# Read-only connection to DuckDB; original Stage 10 models replayed before
# fitting the same endpoints in the canonical boundary-fixed primary cohort.
# Writes only aggregate tables under BPH_WORK_DIR/reporting_rare_stage25.
# Does NOT modify stage outputs, hypotheses, endpoint definitions, or database.

run_stage_25 <- function() {
  options(stringsAsFactors = FALSE, scipen = 999)
  required_packages <- c("DBI", "duckdb", "data.table", "survival", "splines")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) stop("Install packages: ", paste(missing_packages, collapse = ", "))
  suppressPackageStartupMessages({
    library(DBI); library(duckdb); library(data.table)
    library(survival); library(splines)
  })

  R_DIR <- Sys.getenv("BPH_R_DIR", "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R")
  WORK_DIR <- Sys.getenv("BPH_WORK_DIR", "C:/TriNetX/bph_m_work")
  DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")
  OUT_DIR <- file.path(WORK_DIR, "reporting_rare_stage25")
  helper_path <- file.path(R_DIR, "final_analysis_helpers.R")
  stage10_path <- file.path(WORK_DIR, "analysis_final", "10_final_model_results.csv")
  for (f in c(DB_PATH, helper_path, stage10_path)) if (!file.exists(f)) stop("Missing prerequisite: ", f)
  if (dir.exists(OUT_DIR)) stop("Output exists; archive it before rerunning: ", OUT_DIR)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  source(helper_path)
  stage10 <- fread(stage10_path)
  con <- dbConnect(duckdb::duckdb(), dbdir = DB_PATH, read_only = TRUE)
  on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  primary_table <- "analysis_cohort_primary_boundary_fixed"
  if (!primary_table %in% dbListTables(con)) stop("Run Stage 16 first. Missing canonical cohort.")
  current_primary <- as.data.table(dbGetQuery(con,
    "SELECT patient_id FROM analysis_cohort_primary_boundary_fixed"))
  if (nrow(current_primary) != 501445L || anyDuplicated(current_primary$patient_id))
    stop("Canonical cohort is not the verified unique 501,445-patient population")
  setkey(current_primary, patient_id)

  common <- c("patient_id", "exposure", "index_drug", "index_finasteride_strength",
    "index_dutasteride_strength", "index_date", "age_at_index", "race", "ethnicity",
    "marital_status", "patient_regional_location", "patient_source_id", "hypertension",
    "diabetes", "obesity", "ckd", "smoking", "cad", "heart_failure",
    "cerebrovascular_disease", "hyperlipidemia", "osteoporosis", "alopecia",
    "urinary_retention", "uti", "hematuria", "bladder_stone", "n_encounters_prior_365",
    "n_encounters_prior_1095", "n_ambulatory_prior_365", "n_ed_prior_365",
    "n_inpatient_prior_365")
  specs <- data.table(
    source_table = c("analysis_cohort_registry_meningioma",
                     "analysis_cohort_registry_meningioma",
                     "analysis_cohort_treated_meningioma",
                     "analysis_cohort_oncology_treated_registry"),
    original_label = c("ICD-O tumor-table meningioma",
                       "Cancer-registry-source ICD-O meningioma",
                       "Clinically treated meningioma: resection or cranial SRS within 180 d",
                       "Oncology-treatment-table ICD-O meningioma"),
    display_label = c("ICD-O tumor-table meningioma",
                      "Cancer Registry ICD-O meningioma",
                      "Meningioma followed by resection/SRS within 180 d",
                      "Oncology-treatment-table ICD-O meningioma"),
    event_col = c("registry_meningioma_event", "cancer_registry_meningioma_event",
                  "treated_meningioma_event", "oncology_treated_registry_event"),
    time_col = c("registry_followup_days", "cancer_registry_followup_days",
                 "treated_meningioma_followup_days", "oncology_treated_registry_followup_days")
  )

  read_stage10_data <- function(tbl, event, time) {
    fields <- dbListFields(con, tbl)
    mandatory <- c("patient_id", "exposure", "age_at_index", "index_date", event, time)
    if (!all(mandatory %in% fields)) stop(tbl, " missing: ", paste(setdiff(mandatory, fields), collapse = ", "))
    use <- unique(c(common, event, time))
    use <- use[use %in% fields]
    sql <- paste0("SELECT ", paste(sprintf('"%s"', use), collapse = ","), ' FROM "', tbl, '"')
    d <- as.data.table(dbGetQuery(con, sql))
    if (nrow(d) != 501446L || anyDuplicated(d$patient_id))
      stop("Expected one row for each pre-correction primary patient in ", tbl)
    d[, trt := as.integer(exposure == "5ARI")]
    if (anyNA(d$trt) || any(!d$exposure %in% c("5ARI", "URO_ALPHA")))
      stop("Unexpected treatment coding: ", tbl)
    if (anyNA(d[[event]]) || any(!d[[event]] %in% 0:1))
      stop("Invalid event coding in ", tbl, " / ", event)
    if (any(!is.finite(as.numeric(d[[time]]))) || any(as.numeric(d[[time]]) <= 0))
      stop("Nonpositive/missing follow-up in ", tbl, ". Investigate exact-landmark endpoints before refitting.")
    d
  }

  qc <- list(); revised <- list(); revised_incidence <- list();
  for (i in seq_len(nrow(specs))) {
    s <- specs[i]
    cat("Replaying archived model: ", s$original_label, "\n", sep = "")
    if (!s$source_table %in% dbListTables(con)) stop("Missing Stage 09 table: ", s$source_table)
    d <- read_stage10_data(s$source_table, s$event_col, s$time_col)
    archived <- stage10[analysis == s$original_label]
    if (nrow(archived) != 1L) stop("Stage 10 model is absent or duplicated: ", s$original_label)
    replay <- run_ato_generic(d, s$original_label, s$event_col, s$time_col)
    old <- replay$result
    event_ok <- identical(as.integer(old$events[1]), as.integer(archived$events[1]))
    n_ok <- identical(as.integer(old$n[1]), as.integer(archived$n[1]))
    hr_delta <- abs(old$weighted_HR[1] - archived$weighted_HR[1])
    if (!event_ok || !n_ok || !is.finite(hr_delta) || hr_delta > 1e-5)
      stop("Replay does not match archived Stage 10. No revised estimate accepted: ", s$original_label)
    rm(replay); gc(verbose = FALSE)
    setkey(d, patient_id)
    removed <- d[!current_primary, on = "patient_id"]
    if (nrow(removed) != 1L) stop("Rare endpoint table does not differ by exactly one person: ", s$source_table)
    final <- d[current_primary, on = "patient_id", nomatch = 0L]
    if (nrow(final) != 501445L || anyDuplicated(final$patient_id))
      stop("Boundary-aligned rare cohort is not unique N=501445: ", s$original_label)
    result <- run_ato_generic(final, s$original_label, s$event_col, s$time_col)
    rr <- copy(result$result)
    if (rr$n[1] != 501445L) stop("Fitting unexpectedly removed patients: ", s$original_label)
    rr[, `:=`(display_label = s$display_label, source_label = s$original_label)]
    revised[[i]] <- rr
    revised_incidence[[i]] <- result$incidence
    qc[[i]] <- data.table(
      analysis = s$original_label, original_n = old$n[1], revised_n = rr$n[1],
      excluded_n = nrow(removed), original_events = old$events[1],
      revised_events = rr$events[1], original_weighted_HR = old$weighted_HR[1],
      archived_weighted_HR = archived$weighted_HR[1], replay_hr_abs_error = hr_delta,
      revised_weighted_HR = rr$weighted_HR[1], replay_pass = TRUE
    )
    rm(d, final, result); gc(verbose = FALSE)
  }

  models <- rbindlist(revised, fill = TRUE)
  audit <- rbindlist(qc, fill = TRUE)
  incidence <- rbindlist(revised_incidence, fill = TRUE)
  pretty <- models[, .(
    Analysis = display_label,
    N = format(n, big.mark = ",", scientific = FALSE, trim = TRUE),
    Events = events,
    `Adjusted HR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", weighted_HR, weighted_CI_low, weighted_CI_high),
    `P value` = ifelse(weighted_p < 0.001, "<0.001", sprintf("%.3f", weighted_p))
  )]
  fwrite(audit, file.path(OUT_DIR, "25_archived_replay_and_boundary_reconciliation.csv"))
  fwrite(models, file.path(OUT_DIR, "25_rare_boundary_aligned_models_unrounded.csv"))
  fwrite(incidence, file.path(OUT_DIR, "25_rare_boundary_aligned_incidence.csv"))
  fwrite(pretty, file.path(OUT_DIR, "25_Table_S11_Rare_endpoints_BOUNDARY_ALIGNED.csv"))
  writeLines(c(
    "Stage 25 rare endpoint closeout; database read-only; no source tables modified.",
    "Four original Stage 10 model specifications independently replayed and matched before alignment.",
    "Canonical cohort membership adopted by excluding one exact-landmark patient.",
    "Original rare-endpoint definitions and follow-up rules retained; none silently recoded.",
    "Replace Supplementary Table S11 using the generated boundary-aligned CSV only after review.",
    "This operation is a denominator correction, not an independent phenotype revalidation.",
    "Check institutional disclosure policy before sharing aggregate tables externally."
  ), file.path(OUT_DIR, "25_README.txt"))
  message("Stage 25 completed: ", OUT_DIR)
}
run_stage_25()
