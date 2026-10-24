# ============================================================
# 03_refresh_meningioma_outcomes.R
# Finasteride / Dutasteride and Incident Meningioma Study
#
# IMPORTANT CORRECTION:
# Stage 2 used a 365-day landmark but did not explicitly exclude patients
# who developed a meningeal tumor during days 0-364 after treatment initiation.
# Such patients are prevalent cases at the landmark and must not remain at risk.
#
# This script:
#   1) scans diagnosis.csv ONCE for meningeal tumor codes only;
#   2) identifies tumors before index, during lag, and after landmark;
#   3) creates a corrected 1-year-lag primary cohort;
#   4) creates a prespecified 2-year-lag sensitivity cohort;
#   5) creates a >=2-diagnosis confirmation sensitivity endpoint.
#
# Completed heavy work is checkpointed in DuckDB and reused on rerun.
# ============================================================

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
PARQUET_DIR <- file.path(WORK_DIR, "parquet")
DIAGNOSIS_FILE <- file.path(RAW_DIR, "diagnosis.csv")

dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PARQUET_DIR, recursive = TRUE, showWarnings = FALSE)

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

if (!"analysis_cohort_primary" %in% dbListTables(con)) {
  stop("analysis_cohort_primary not found. Run Stage 2 v3 first.")
}

# ------------------------------------------------------------
# 1. MENINGIOMA-ONLY RAW SCAN
# ------------------------------------------------------------
# Strict primary phenotype:
#   ICD-10-CM D32.0 intracranial meninges
#   ICD-10-CM D32.9 meninges, unspecified
#   ICD-9-CM 225.2  benign neoplasm of cerebral meninges
#
# Broad exclusion/sensitivity phenotype:
#   strict +
#   D42.0/D42.9 uncertain behavior meninges
#   C70.0/C70.9 malignant meninges
#   ICD-9 237.6 / 192.1
#
# We use the BROAD phenotype during the lag exclusion so a suspicious
# meningeal-tumor code during the washout cannot later masquerade as incident
# benign meningioma.
if (REUSE_EXISTING_HEAVY_TABLES &&
    "meningioma_dx_events" %in% dbListTables(con)) {

  cat("\nReusing existing meningioma_dx_events; diagnosis.csv will NOT be rescanned.\n")

} else {

  cat("\n============================================================\n")
  cat("Scanning diagnosis.csv for meningeal tumor codes only...\n")
  cat("============================================================\n")
  t0 <- Sys.time()

  dx_sql <- paste0(
"
CREATE OR REPLACE TABLE meningioma_dx_events AS
WITH d AS (
    SELECT
        patient_id,
        UPPER(TRIM(code_system)) AS code_system,
        UPPER(REPLACE(TRIM(code), '.', '')) AS code_norm,
        TRY_CAST(STRPTIME(date, '%Y%m%d') AS DATE) AS dx_date
    FROM read_csv_auto(
        '", sql_path(DIAGNOSIS_FILE), "',
        header=true,
        all_varchar=true,
        ignore_errors=false
    )
    WHERE date IS NOT NULL
      AND date <> ''
      AND code IS NOT NULL
      AND code <> ''
),
m AS (
    SELECT
        patient_id,
        dx_date,
        code_system,
        code_norm,

        CASE
          WHEN code_system='ICD-10-CM'
           AND code_norm IN ('D320','D329')
            THEN 1
          WHEN code_system='ICD-9-CM'
           AND code_norm='2252'
            THEN 1
          ELSE 0
        END AS strict_meningioma,

        CASE
          WHEN code_system='ICD-10-CM'
           AND code_norm IN ('D320','D329','D420','D429','C700','C709')
            THEN 1
          WHEN code_system='ICD-9-CM'
           AND code_norm IN ('2252','2376','1921')
            THEN 1
          ELSE 0
        END AS broad_meningeal_tumor

    FROM d
    WHERE
      (
        code_system='ICD-10-CM'
        AND code_norm IN ('D320','D329','D420','D429','C700','C709')
      )
      OR
      (
        code_system='ICD-9-CM'
        AND code_norm IN ('2252','2376','1921')
      )
)
SELECT DISTINCT
    m.patient_id,
    m.dx_date,
    m.code_system,
    m.code_norm,
    m.strict_meningioma,
    m.broad_meningeal_tumor
FROM m
INNER JOIN analysis_cohort_primary a
  ON a.patient_id = m.patient_id
WHERE m.dx_date IS NOT NULL;
"
  )

  dbExecute(con, dx_sql)

  try(dbExecute(
    con,
    "CREATE INDEX IF NOT EXISTS idx_meningioma_events_patient_date
     ON meningioma_dx_events(patient_id, dx_date);"
  ), silent = TRUE)

  cat(
    "Meningioma-only diagnosis scan completed in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    " minutes.\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# 2. PATIENT-LEVEL TEMPORAL OUTCOME SUMMARY
# ------------------------------------------------------------
dbExecute(con, "
CREATE OR REPLACE TABLE meningioma_temporal_summary AS
WITH x AS (
    SELECT
        a.patient_id,
        a.index_date,
        a.last_encounter_date,
        m.dx_date,
        m.strict_meningioma,
        m.broad_meningeal_tumor
    FROM analysis_cohort_primary a
    LEFT JOIN meningioma_dx_events m
      ON m.patient_id = a.patient_id
),
base AS (
    SELECT
        patient_id,

        MIN(dx_date) FILTER (
          WHERE broad_meningeal_tumor=1
            AND dx_date < index_date
        ) AS broad_preindex_date,

        MIN(dx_date) FILTER (
          WHERE strict_meningioma=1
            AND dx_date < index_date
        ) AS strict_preindex_date,

        MIN(dx_date) FILTER (
          WHERE broad_meningeal_tumor=1
            AND dx_date >= index_date
            AND dx_date < index_date + INTERVAL 365 DAY
        ) AS broad_lag1_date,

        MIN(dx_date) FILTER (
          WHERE strict_meningioma=1
            AND dx_date >= index_date
            AND dx_date < index_date + INTERVAL 365 DAY
        ) AS strict_lag1_date,

        MIN(dx_date) FILTER (
          WHERE broad_meningeal_tumor=1
            AND dx_date >= index_date
            AND dx_date < index_date + INTERVAL 730 DAY
        ) AS broad_lag2_date,

        MIN(dx_date) FILTER (
          WHERE strict_meningioma=1
            AND dx_date >= index_date
            AND dx_date < index_date + INTERVAL 730 DAY
        ) AS strict_lag2_date,

        MIN(dx_date) FILTER (
          WHERE strict_meningioma=1
            AND dx_date >= index_date + INTERVAL 365 DAY
            AND dx_date <= last_encounter_date
        ) AS strict_post_lag1_first,

        MIN(dx_date) FILTER (
          WHERE broad_meningeal_tumor=1
            AND dx_date >= index_date + INTERVAL 365 DAY
            AND dx_date <= last_encounter_date
        ) AS broad_post_lag1_first,

        MIN(dx_date) FILTER (
          WHERE strict_meningioma=1
            AND dx_date >= index_date + INTERVAL 730 DAY
            AND dx_date <= last_encounter_date
        ) AS strict_post_lag2_first,

        MIN(dx_date) FILTER (
          WHERE broad_meningeal_tumor=1
            AND dx_date >= index_date + INTERVAL 730 DAY
            AND dx_date <= last_encounter_date
        ) AS broad_post_lag2_first,

        COUNT(DISTINCT dx_date) FILTER (
          WHERE strict_meningioma=1
            AND dx_date >= index_date + INTERVAL 365 DAY
            AND dx_date <= last_encounter_date
        ) AS n_strict_dates_post_lag1

    FROM x
    GROUP BY patient_id
),
confirm1 AS (
    SELECT
        b.patient_id,
        MIN(m.dx_date) AS strict_lag1_confirmation_date
    FROM base b
    JOIN meningioma_dx_events m
      ON m.patient_id=b.patient_id
     AND m.strict_meningioma=1
     AND b.strict_post_lag1_first IS NOT NULL
     AND m.dx_date >= b.strict_post_lag1_first + INTERVAL 30 DAY
    JOIN analysis_cohort_primary a
      ON a.patient_id=b.patient_id
     AND m.dx_date <= a.last_encounter_date
    GROUP BY b.patient_id
),
confirm2 AS (
    SELECT
        b.patient_id,
        MIN(m.dx_date) AS strict_lag2_confirmation_date
    FROM base b
    JOIN meningioma_dx_events m
      ON m.patient_id=b.patient_id
     AND m.strict_meningioma=1
     AND b.strict_post_lag2_first IS NOT NULL
     AND m.dx_date >= b.strict_post_lag2_first + INTERVAL 30 DAY
    JOIN analysis_cohort_primary a
      ON a.patient_id=b.patient_id
     AND m.dx_date <= a.last_encounter_date
    GROUP BY b.patient_id
)
SELECT
    b.*,
    c1.strict_lag1_confirmation_date,
    c2.strict_lag2_confirmation_date
FROM base b
LEFT JOIN confirm1 c1 USING (patient_id)
LEFT JOIN confirm2 c2 USING (patient_id);
")

# ------------------------------------------------------------
# 3. CORRECTED 1-YEAR LANDMARK COHORT
# ------------------------------------------------------------
# Exclude ANY broad meningeal tumor during the 365-day lag.
dbExecute(con, "
CREATE OR REPLACE TABLE analysis_cohort_primary_corrected AS
SELECT
    a.* EXCLUDE (
      incident_meningioma_date_strict,
      incident_meningeal_tumor_date_broad,
      meningioma_event_strict,
      meningeal_event_broad,
      followup_end_strict,
      followup_end_broad,
      followup_days_strict,
      followup_days_broad
    ),

    t.strict_post_lag1_first AS incident_meningioma_date_strict,
    t.broad_post_lag1_first AS incident_meningeal_tumor_date_broad,

    CASE
      WHEN t.strict_post_lag1_first IS NOT NULL THEN 1 ELSE 0
    END AS meningioma_event_strict,

    CASE
      WHEN t.broad_post_lag1_first IS NOT NULL THEN 1 ELSE 0
    END AS meningeal_event_broad,

    CASE
      WHEN t.strict_lag1_confirmation_date IS NOT NULL THEN 1 ELSE 0
    END AS meningioma_event_confirmed,

    t.strict_lag1_confirmation_date AS meningioma_confirmation_date,

    CASE
      WHEN t.strict_post_lag1_first IS NOT NULL
      THEN t.strict_post_lag1_first
      ELSE a.last_encounter_date
    END AS followup_end_strict,

    DATE_DIFF(
      'day',
      a.index_date + INTERVAL 365 DAY,
      CASE
        WHEN t.strict_post_lag1_first IS NOT NULL
        THEN t.strict_post_lag1_first
        ELSE a.last_encounter_date
      END
    ) AS followup_days_strict

FROM analysis_cohort_primary a
JOIN meningioma_temporal_summary t USING (patient_id)
WHERE t.broad_preindex_date IS NULL
  AND t.broad_lag1_date IS NULL;
")

try(dbExecute(
  con,
  "CREATE INDEX IF NOT EXISTS idx_primary_corrected_patient
   ON analysis_cohort_primary_corrected(patient_id);"
), silent = TRUE)

# ------------------------------------------------------------
# 4. 2-YEAR LAG SENSITIVITY COHORT
# ------------------------------------------------------------
# Need >=365 days of observation after the 2-year landmark.
dbExecute(con, "
CREATE OR REPLACE TABLE analysis_cohort_lag2 AS
SELECT
    a.* EXCLUDE (
      incident_meningioma_date_strict,
      incident_meningeal_tumor_date_broad,
      meningioma_event_strict,
      meningeal_event_broad,
      followup_end_strict,
      followup_end_broad,
      followup_days_strict,
      followup_days_broad
    ),

    a.index_date + INTERVAL 730 DAY AS landmark_date_lag2,

    t.strict_post_lag2_first AS incident_meningioma_date_strict,

    CASE
      WHEN t.strict_post_lag2_first IS NOT NULL THEN 1 ELSE 0
    END AS meningioma_event_strict,

    CASE
      WHEN t.strict_lag2_confirmation_date IS NOT NULL THEN 1 ELSE 0
    END AS meningioma_event_confirmed,

    t.strict_lag2_confirmation_date AS meningioma_confirmation_date,

    CASE
      WHEN t.strict_post_lag2_first IS NOT NULL
      THEN t.strict_post_lag2_first
      ELSE a.last_encounter_date
    END AS followup_end_strict,

    DATE_DIFF(
      'day',
      a.index_date + INTERVAL 730 DAY,
      CASE
        WHEN t.strict_post_lag2_first IS NOT NULL
        THEN t.strict_post_lag2_first
        ELSE a.last_encounter_date
      END
    ) AS followup_days_strict

FROM analysis_cohort_primary a
JOIN meningioma_temporal_summary t USING (patient_id)
WHERE t.broad_preindex_date IS NULL
  AND t.broad_lag2_date IS NULL
  AND a.last_encounter_date >= a.index_date + INTERVAL 1095 DAY;
")

# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------
lag_exclusions <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    COUNT(*) AS stage2_patients,
    SUM(CASE WHEN t.broad_lag1_date IS NOT NULL THEN 1 ELSE 0 END) AS broad_tumor_during_first_365d,
    SUM(CASE WHEN t.strict_lag1_date IS NOT NULL THEN 1 ELSE 0 END) AS strict_meningioma_during_first_365d,
    SUM(CASE WHEN t.broad_lag2_date IS NOT NULL THEN 1 ELSE 0 END) AS broad_tumor_during_first_730d,
    SUM(CASE WHEN t.strict_lag2_date IS NOT NULL THEN 1 ELSE 0 END) AS strict_meningioma_during_first_730d
FROM analysis_cohort_primary a
JOIN meningioma_temporal_summary t USING (patient_id)
GROUP BY exposure
ORDER BY exposure;
"))

primary_qc <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    COUNT(*) AS patients,
    SUM(meningioma_event_strict) AS strict_events,
    ROUND(100.0*SUM(meningioma_event_strict)/COUNT(*),4) AS strict_event_pct,
    SUM(meningioma_event_confirmed) AS confirmed_events,
    ROUND(100.0*SUM(meningioma_event_confirmed)/COUNT(*),4) AS confirmed_event_pct,
    SUM(followup_days_strict)/365.25 AS total_person_years,
    100000.0*SUM(meningioma_event_strict)/(SUM(followup_days_strict)/365.25)
      AS crude_incidence_per_100k_py,
    MEDIAN(followup_days_strict)/365.25 AS median_followup_y
FROM analysis_cohort_primary_corrected
GROUP BY exposure
ORDER BY exposure;
"))

lag2_qc <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    COUNT(*) AS patients,
    SUM(meningioma_event_strict) AS strict_events,
    SUM(meningioma_event_confirmed) AS confirmed_events,
    SUM(followup_days_strict)/365.25 AS total_person_years,
    100000.0*SUM(meningioma_event_strict)/(SUM(followup_days_strict)/365.25)
      AS crude_incidence_per_100k_py,
    MEDIAN(followup_days_strict)/365.25 AS median_followup_y
FROM analysis_cohort_lag2
GROUP BY exposure
ORDER BY exposure;
"))

fwrite(lag_exclusions, file.path(QC_DIR, "03_lag_exclusions.csv"))
fwrite(primary_qc, file.path(QC_DIR, "03_primary_corrected_counts.csv"))
fwrite(lag2_qc, file.path(QC_DIR, "03_lag2_counts.csv"))

# ------------------------------------------------------------
# 6. EXPORT
# ------------------------------------------------------------
p1 <- sql_path(file.path(PARQUET_DIR, "analysis_cohort_primary_corrected.parquet"))
p2 <- sql_path(file.path(PARQUET_DIR, "analysis_cohort_lag2.parquet"))

dbExecute(con, sprintf(
  "COPY analysis_cohort_primary_corrected TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000);", p1
))
dbExecute(con, sprintf(
  "COPY analysis_cohort_lag2 TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000);", p2
))

cat("\n================ LAG EXCLUSIONS ====================\n")
print(lag_exclusions)

cat("\n================ CORRECTED PRIMARY =================\n")
print(primary_qc)

cat("\n================ 2-YEAR LAG ========================\n")
print(lag2_qc)

cat("\nStage 3 complete.\n")

dbDisconnect(con, shutdown = TRUE)
