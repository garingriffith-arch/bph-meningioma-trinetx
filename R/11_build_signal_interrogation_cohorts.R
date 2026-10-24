# ============================================================
# 11_build_signal_interrogation_cohorts.R
#
# Targeted follow-up to the positive exploratory findings:
#   A) sustained/persistent exposure + healthcare utilization during
#      the exposure-assessment landmark period;
#   B) finasteride- and dutasteride-specific sustained cohorts;
#   C) contemporary US-only exact-dose finasteride vs dutasteride cohorts.
#
# New heavy scan:
#   encounter.csv ONCE, checkpointed.
#
# IMPORTANT:
# Post-index healthcare utilization is measured ONLY during the
# exposure-assessment period and outcome follow-up begins after that period.
# These models are bias-probing sensitivities, not replacements for the
# primary new-user analysis.
# ============================================================

run_stage_11 <- function() {

  options(stringsAsFactors = FALSE, scipen = 999)

  required_pkgs <- c("DBI", "duckdb", "data.table")
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

  library(DBI)
  library(duckdb)
  library(data.table)

  RAW_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/Data/BPH"
  WORK_DIR <- "C:/TriNetX/bph_m_work"
  DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")
  TEMP_DIR <- file.path(WORK_DIR, "duckdb_temp")
  QC_DIR <- file.path(WORK_DIR, "qc")
  ENCOUNTER_FILE <- file.path(RAW_DIR, "encounter.csv")

  dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

  N_THREADS <- min(24L, max(4L, parallel::detectCores(logical = TRUE) - 6L))
  MEMORY_LIMIT <- "48GB"
  REUSE_EXISTING_HEAVY_TABLES <- TRUE

  sql_path <- function(x) {
    x <- normalizePath(x, winslash = "/", mustWork = FALSE)
    gsub("'", "''", x, fixed = TRUE)
  }

  con <- dbConnect(duckdb::duckdb(), dbdir = DB_PATH, read_only = FALSE)
  on.exit({
    try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }, add = TRUE)

  dbExecute(con, sprintf("PRAGMA threads=%d;", N_THREADS))
  dbExecute(con, sprintf("PRAGMA memory_limit='%s';", MEMORY_LIMIT))
  dbExecute(con, sprintf("PRAGMA temp_directory='%s';", sql_path(TEMP_DIR)))
  dbExecute(con, "SET preserve_insertion_order=false;")

  needed <- c(
    "analysis_cohort_primary_revised",
    "analysis_cohort_sustained_1y",
    "analysis_cohort_sustained_1y_strict",
    "analysis_cohort_persistent_3y"
  )
  missing <- setdiff(needed, dbListTables(con))
  if (length(missing) > 0L) {
    stop(
      "Run the final sensitivity suite through Stage 10 first. Missing: ",
      paste(missing, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 1. TARGET PATIENTS FOR LANDMARK UTILIZATION
  # ----------------------------------------------------------
  dbExecute(con, "
  CREATE OR REPLACE TABLE signal_interrogation_targets AS
  SELECT
    patient_id,
    index_date,
    365 AS landmark_days
  FROM analysis_cohort_sustained_1y

  UNION

  SELECT
    patient_id,
    index_date,
    1095 AS landmark_days
  FROM analysis_cohort_persistent_3y;
  ")

  # ----------------------------------------------------------
  # 2. ONE ENCOUNTER SCAN
  # ----------------------------------------------------------
  if (
    REUSE_EXISTING_HEAVY_TABLES &&
    "signal_landmark_utilization" %in% dbListTables(con)
  ) {

    cat(
      "\nReusing signal_landmark_utilization; ",
      "encounter.csv will NOT be rescanned.\n",
      sep = ""
    )

  } else {

    cat("\n============================================================\n")
    cat("Scanning encounter.csv for landmark-period utilization...\n")
    cat("============================================================\n")
    t0 <- Sys.time()

    enc_sql <- paste0(
"
CREATE OR REPLACE TABLE signal_landmark_utilization AS
WITH e AS (
  SELECT
    patient_id,
    TRY_CAST(STRPTIME(start_date, '%Y%m%d') AS DATE) AS encounter_date,
    UPPER(TRIM(type)) AS encounter_type
  FROM read_csv_auto(
    '", sql_path(ENCOUNTER_FILE), "',
    header=true,
    all_varchar=true,
    ignore_errors=false
  )
  WHERE start_date IS NOT NULL
    AND start_date <> ''
),
j AS (
  SELECT
    t.patient_id,
    t.index_date,
    t.landmark_days,
    e.encounter_date,
    e.encounter_type
  FROM signal_interrogation_targets t
  LEFT JOIN e
    ON e.patient_id=t.patient_id
   AND e.encounter_date >= t.index_date
   AND e.encounter_date < t.index_date + t.landmark_days * INTERVAL 1 DAY
)
SELECT
  patient_id,
  landmark_days,

  COUNT(encounter_date) AS lm_encounters_total,
  COUNT(DISTINCT encounter_date) AS lm_encounter_days,

  COUNT(encounter_date) FILTER (
    WHERE encounter_type='AMB'
  ) AS lm_ambulatory,

  COUNT(encounter_date) FILTER (
    WHERE encounter_type='EMER'
  ) AS lm_ed,

  COUNT(encounter_date) FILTER (
    WHERE encounter_type IN ('IMP','ACUTE','NONAC','SS')
  ) AS lm_inpatient,

  COUNT(DISTINCT encounter_date) FILTER (
    WHERE encounter_date >= index_date
      AND encounter_date < index_date + INTERVAL 365 DAY
  ) AS lm_y1_encounter_days,

  COUNT(DISTINCT encounter_date) FILTER (
    WHERE landmark_days >= 730
      AND encounter_date >= index_date + INTERVAL 365 DAY
      AND encounter_date < index_date + INTERVAL 730 DAY
  ) AS lm_y2_encounter_days,

  COUNT(DISTINCT encounter_date) FILTER (
    WHERE landmark_days >= 1095
      AND encounter_date >= index_date + INTERVAL 730 DAY
      AND encounter_date < index_date + INTERVAL 1095 DAY
  ) AS lm_y3_encounter_days

FROM j
GROUP BY patient_id, landmark_days;
"
    )

    dbExecute(con, enc_sql)

    try(
      dbExecute(
        con,
        "CREATE INDEX IF NOT EXISTS idx_signal_lm_util
         ON signal_landmark_utilization(patient_id, landmark_days);"
      ),
      silent = TRUE
    )

    cat(
      "Landmark utilization scan completed in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      " minutes.\n",
      sep = ""
    )
  }

  # ----------------------------------------------------------
  # 3. ENRICH SUSTAINED / PERSISTENT COHORTS
  # ----------------------------------------------------------
  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_sustained_1y_util AS
  SELECT
    a.*,
    COALESCE(u.lm_encounters_total,0) AS lm_encounters_total,
    COALESCE(u.lm_encounter_days,0) AS lm_encounter_days,
    COALESCE(u.lm_ambulatory,0) AS lm_ambulatory,
    COALESCE(u.lm_ed,0) AS lm_ed,
    COALESCE(u.lm_inpatient,0) AS lm_inpatient,
    COALESCE(u.lm_y1_encounter_days,0) AS lm_y1_encounter_days
  FROM analysis_cohort_sustained_1y a
  LEFT JOIN signal_landmark_utilization u
    ON u.patient_id=a.patient_id
   AND u.landmark_days=365;
  ")

  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_sustained_1y_strict_util AS
  SELECT
    a.*,
    COALESCE(u.lm_encounters_total,0) AS lm_encounters_total,
    COALESCE(u.lm_encounter_days,0) AS lm_encounter_days,
    COALESCE(u.lm_ambulatory,0) AS lm_ambulatory,
    COALESCE(u.lm_ed,0) AS lm_ed,
    COALESCE(u.lm_inpatient,0) AS lm_inpatient,
    COALESCE(u.lm_y1_encounter_days,0) AS lm_y1_encounter_days
  FROM analysis_cohort_sustained_1y_strict a
  LEFT JOIN signal_landmark_utilization u
    ON u.patient_id=a.patient_id
   AND u.landmark_days=365;
  ")

  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_persistent_3y_util AS
  SELECT
    a.*,
    COALESCE(u.lm_encounters_total,0) AS lm_encounters_total,
    COALESCE(u.lm_encounter_days,0) AS lm_encounter_days,
    COALESCE(u.lm_ambulatory,0) AS lm_ambulatory,
    COALESCE(u.lm_ed,0) AS lm_ed,
    COALESCE(u.lm_inpatient,0) AS lm_inpatient,
    COALESCE(u.lm_y1_encounter_days,0) AS lm_y1_encounter_days,
    COALESCE(u.lm_y2_encounter_days,0) AS lm_y2_encounter_days,
    COALESCE(u.lm_y3_encounter_days,0) AS lm_y3_encounter_days
  FROM analysis_cohort_persistent_3y a
  LEFT JOIN signal_landmark_utilization u
    ON u.patient_id=a.patient_id
   AND u.landmark_days=1095;
  ")

  # ----------------------------------------------------------
  # 4. REGION AUDIT
  # ----------------------------------------------------------
  region_qc <- as.data.table(
    dbGetQuery(
      con,
      "
      SELECT
        patient_regional_location,
        index_drug,
        COUNT(*) AS patients,
        SUM(meningioma_event_strict) AS events
      FROM analysis_cohort_primary_revised
      WHERE exposure='5ARI'
        AND index_drug IN ('finasteride','dutasteride')
      GROUP BY patient_regional_location, index_drug
      ORDER BY patient_regional_location, index_drug;
      "
    )
  )

  fwrite(
    region_qc,
    file.path(QC_DIR, "11_head_to_head_region_audit.csv")
  )

  # TriNetX regional coding in this export contains "Ex-US", "Unknown", and
  # US Census divisions. Therefore a defensible US-only restriction is:
  # observed region != Ex-US and != Unknown.
  #
  # Do not silently include missing geography.
  us_condition <- "
    patient_regional_location IS NOT NULL
    AND TRIM(patient_regional_location) <> ''
    AND UPPER(TRIM(patient_regional_location)) NOT IN ('EX-US','UNKNOWN')
  "

  # ----------------------------------------------------------
  # 5. CONTEMPORARY US EXACT-DOSE HEAD-TO-HEAD COHORTS
  # ----------------------------------------------------------
  dbExecute(
    con,
    paste0(
"
CREATE OR REPLACE TABLE analysis_cohort_fd_us_2012_exact AS
SELECT *
FROM analysis_cohort_primary_revised
WHERE exposure='5ARI'
  AND index_date >= DATE '2012-01-01'
  AND ", us_condition, "
  AND (
    (
      index_drug='finasteride'
      AND UPPER(TRIM(COALESCE(index_finasteride_strength,'')))='5 MG'
    )
    OR
    (
      index_drug='dutasteride'
      AND UPPER(TRIM(COALESCE(index_dutasteride_strength,'')))='0.5 MG'
    )
  );
"
    )
  )

  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_fd_us_2012_exact_age50 AS
  SELECT *
  FROM analysis_cohort_fd_us_2012_exact
  WHERE age_at_index >= 50;
  ")

  dbExecute(
    con,
    paste0(
"
CREATE OR REPLACE TABLE analysis_cohort_fd_us_icd10_exact AS
SELECT *
FROM analysis_cohort_primary_revised
WHERE exposure='5ARI'
  AND index_date >= DATE '2015-10-01'
  AND ", us_condition, "
  AND (
    (
      index_drug='finasteride'
      AND UPPER(TRIM(COALESCE(index_finasteride_strength,'')))='5 MG'
    )
    OR
    (
      index_drug='dutasteride'
      AND UPPER(TRIM(COALESCE(index_dutasteride_strength,'')))='0.5 MG'
    )
  );
"
    )
  )

  # ----------------------------------------------------------
  # 6. QC TABLES
  # ----------------------------------------------------------
  util_qc <- as.data.table(
    dbGetQuery(
      con,
      "
      SELECT
        'sustained_1y' AS cohort,
        exposure,
        COUNT(*) AS patients,
        SUM(meningioma_event_strict) AS events,
        MEDIAN(lm_encounter_days) AS median_landmark_encounter_days,
        MEDIAN(lm_ambulatory) AS median_landmark_ambulatory,
        MEDIAN(lm_ed) AS median_landmark_ed,
        MEDIAN(lm_inpatient) AS median_landmark_inpatient
      FROM analysis_cohort_sustained_1y_util
      GROUP BY exposure

      UNION ALL

      SELECT
        'sustained_1y_strict' AS cohort,
        exposure,
        COUNT(*) AS patients,
        SUM(meningioma_event_strict) AS events,
        MEDIAN(lm_encounter_days) AS median_landmark_encounter_days,
        MEDIAN(lm_ambulatory) AS median_landmark_ambulatory,
        MEDIAN(lm_ed) AS median_landmark_ed,
        MEDIAN(lm_inpatient) AS median_landmark_inpatient
      FROM analysis_cohort_sustained_1y_strict_util
      GROUP BY exposure

      UNION ALL

      SELECT
        'persistent_3y' AS cohort,
        exposure,
        COUNT(*) AS patients,
        SUM(meningioma_event_strict) AS events,
        MEDIAN(lm_encounter_days) AS median_landmark_encounter_days,
        MEDIAN(lm_ambulatory) AS median_landmark_ambulatory,
        MEDIAN(lm_ed) AS median_landmark_ed,
        MEDIAN(lm_inpatient) AS median_landmark_inpatient
      FROM analysis_cohort_persistent_3y_util
      GROUP BY exposure

      ORDER BY cohort, exposure;
      "
    )
  )

  fd_qc <- as.data.table(
    dbGetQuery(
      con,
      "
      SELECT
        'US_2012_exact' AS cohort,
        index_drug,
        COUNT(*) AS patients,
        SUM(meningioma_event_strict) AS events
      FROM analysis_cohort_fd_us_2012_exact
      GROUP BY index_drug

      UNION ALL

      SELECT
        'US_2012_exact_age50' AS cohort,
        index_drug,
        COUNT(*) AS patients,
        SUM(meningioma_event_strict) AS events
      FROM analysis_cohort_fd_us_2012_exact_age50
      GROUP BY index_drug

      UNION ALL

      SELECT
        'US_ICD10_exact' AS cohort,
        index_drug,
        COUNT(*) AS patients,
        SUM(meningioma_event_strict) AS events
      FROM analysis_cohort_fd_us_icd10_exact
      GROUP BY index_drug

      ORDER BY cohort, index_drug;
      "
    )
  )

  fwrite(
    util_qc,
    file.path(QC_DIR, "11_landmark_utilization_qc.csv")
  )
  fwrite(
    fd_qc,
    file.path(QC_DIR, "11_fd_homogeneous_cohort_counts.csv")
  )

  cat("\n================ LANDMARK UTILIZATION QC ============\n")
  print(util_qc)

  cat("\n================ FINASTERIDE VS DUTASTERIDE QC ======\n")
  print(fd_qc)

  cat("\nStage 11 complete.\n")
}

run_stage_11()
