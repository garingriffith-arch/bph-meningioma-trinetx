# ============================================================
# 05_rebuild_without_future_followup_requirement.R
# Finasteride / Dutasteride and Incident Meningioma Study
#
# DESIGN CORRECTION
# -----------------
# The prior Stage 2 primary cohort required:
#   last encounter >= index + 730 days
#
# Because risk follow-up begins at index + 365 days, that criterion selects
# patients based on a FULL YEAR of future observed healthcare after the
# landmark. It excluded >100,000 otherwise eligible patients and conditions
# inclusion on future utilization.
#
# Revised primary design:
#   - keep all prior eligibility criteria;
#   - require only that the patient has some observed follow-up AFTER the
#     365-day landmark;
#   - censor at the patient's last observed encounter;
#   - preserve the original index calendar window for comparability;
#   - exclude any broad meningeal tumor before or during the 365-day lag.
#
# A revised 2-year-lag sensitivity cohort is built analogously.
#
# Heavy scan:
#   diagnosis.csv once, but only meningeal-tumor codes are retained.
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

if (!"cohort_eligibility" %in% dbListTables(con)) {
  stop("cohort_eligibility not found. Run Stage 2 v3 first.")
}

# ------------------------------------------------------------
# 1. REVISED ELIGIBLE BASE
# ------------------------------------------------------------
# Preserve every Stage-2 restriction EXCEPT the requirement for a full
# additional year of future observation after the landmark.
#
# Require last encounter strictly after the landmark so each patient has
# positive observable risk time.
dbExecute(con, "
CREATE OR REPLACE TABLE revised_eligible_base AS
SELECT *
FROM cohort_eligibility
WHERE eligible_male = 1
  AND eligible_age = 1
  AND eligible_calendar = 1
  AND eligible_lookback = 1
  AND eligible_bph = 1
  AND eligible_no_prostate_cancer = 1
  AND eligible_no_prior_meningeal_tumor = 1
  AND opposite_class_within_30d = 0
  AND eligible_not_known_finasteride_1mg = 1
  AND last_encounter_date > index_date + INTERVAL 365 DAY;
")

base_counts <- as.data.table(dbGetQuery(con, "
SELECT
  exposure,
  COUNT(*) AS patients,
  MEDIAN(DATE_DIFF('day', index_date + INTERVAL 365 DAY, last_encounter_date))/365.25
    AS median_potential_post_landmark_followup_y
FROM revised_eligible_base
GROUP BY exposure
ORDER BY exposure;
"))

fwrite(base_counts, file.path(QC_DIR, "05_revised_base_counts.csv"))

# ------------------------------------------------------------
# 2. MENINGIOMA-ONLY DIAGNOSIS SCAN
# ------------------------------------------------------------
if (REUSE_EXISTING_HEAVY_TABLES &&
    "meningioma_dx_events_revised" %in% dbListTables(con)) {

  cat("\nReusing meningioma_dx_events_revised; diagnosis.csv will NOT be rescanned.\n")

} else {

  cat("\n============================================================\n")
  cat("Scanning diagnosis.csv for revised meningioma cohort...\n")
  cat("============================================================\n")
  t0 <- Sys.time()

  dx_sql <- paste0(
"
CREATE OR REPLACE TABLE meningioma_dx_events_revised AS
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
INNER JOIN revised_eligible_base b
  ON b.patient_id = m.patient_id
WHERE m.dx_date IS NOT NULL;
"
  )

  dbExecute(con, dx_sql)

  try(dbExecute(
    con,
    "CREATE INDEX IF NOT EXISTS idx_meningioma_events_revised_patient_date
     ON meningioma_dx_events_revised(patient_id, dx_date);"
  ), silent = TRUE)

  cat(
    "Revised meningioma diagnosis scan completed in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    " minutes.\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# 3. TEMPORAL SUMMARY
# ------------------------------------------------------------
dbExecute(con, "
CREATE OR REPLACE TABLE meningioma_temporal_summary_revised AS
WITH x AS (
    SELECT
        a.patient_id,
        a.index_date,
        a.last_encounter_date,
        m.dx_date,
        m.strict_meningioma,
        m.broad_meningeal_tumor
    FROM revised_eligible_base a
    LEFT JOIN meningioma_dx_events_revised m
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
          WHERE broad_meningeal_tumor=1
            AND dx_date >= index_date
            AND dx_date < index_date + INTERVAL 365 DAY
        ) AS broad_lag1_date,

        MIN(dx_date) FILTER (
          WHERE broad_meningeal_tumor=1
            AND dx_date >= index_date
            AND dx_date < index_date + INTERVAL 730 DAY
        ) AS broad_lag2_date,

        MIN(dx_date) FILTER (
          WHERE strict_meningioma=1
            AND dx_date >= index_date + INTERVAL 365 DAY
            AND dx_date <= last_encounter_date
        ) AS strict_post_lag1_first,

        MIN(dx_date) FILTER (
          WHERE strict_meningioma=1
            AND dx_date >= index_date + INTERVAL 730 DAY
            AND dx_date <= last_encounter_date
        ) AS strict_post_lag2_first,

        MIN(dx_date) FILTER (
          WHERE broad_meningeal_tumor=1
            AND dx_date >= index_date + INTERVAL 365 DAY
            AND dx_date <= last_encounter_date
        ) AS broad_post_lag1_first,

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
    JOIN meningioma_dx_events_revised m
      ON m.patient_id=b.patient_id
     AND m.strict_meningioma=1
     AND b.strict_post_lag1_first IS NOT NULL
     AND m.dx_date >= b.strict_post_lag1_first + INTERVAL 30 DAY
    JOIN revised_eligible_base a
      ON a.patient_id=b.patient_id
     AND m.dx_date <= a.last_encounter_date
    GROUP BY b.patient_id
),
confirm2 AS (
    SELECT
        b.patient_id,
        MIN(m.dx_date) AS strict_lag2_confirmation_date
    FROM base b
    JOIN meningioma_dx_events_revised m
      ON m.patient_id=b.patient_id
     AND m.strict_meningioma=1
     AND b.strict_post_lag2_first IS NOT NULL
     AND m.dx_date >= b.strict_post_lag2_first + INTERVAL 30 DAY
    JOIN revised_eligible_base a
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
# 4. REVISED PRIMARY COHORT
# ------------------------------------------------------------
dbExecute(con, "
CREATE OR REPLACE TABLE analysis_cohort_primary_revised AS
SELECT
    a.*,

    a.index_date + INTERVAL 365 DAY AS landmark_date,

    t.strict_post_lag1_first AS incident_meningioma_date_strict,
    t.broad_post_lag1_first AS incident_meningeal_tumor_date_broad,

    CASE WHEN t.strict_post_lag1_first IS NOT NULL THEN 1 ELSE 0 END
      AS meningioma_event_strict,

    CASE WHEN t.broad_post_lag1_first IS NOT NULL THEN 1 ELSE 0 END
      AS meningeal_event_broad,

    CASE WHEN t.strict_lag1_confirmation_date IS NOT NULL THEN 1 ELSE 0 END
      AS meningioma_event_confirmed,

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

FROM revised_eligible_base a
JOIN meningioma_temporal_summary_revised t USING (patient_id)
WHERE t.broad_preindex_date IS NULL
  AND t.broad_lag1_date IS NULL
  AND a.last_encounter_date > a.index_date + INTERVAL 365 DAY;
")

# ------------------------------------------------------------
# 5. REVISED 2-YEAR-LAG SENSITIVITY COHORT
# ------------------------------------------------------------
dbExecute(con, "
CREATE OR REPLACE TABLE analysis_cohort_lag2_revised AS
SELECT
    a.*,

    a.index_date + INTERVAL 730 DAY AS landmark_date,

    t.strict_post_lag2_first AS incident_meningioma_date_strict,

    CASE WHEN t.strict_post_lag2_first IS NOT NULL THEN 1 ELSE 0 END
      AS meningioma_event_strict,

    CASE WHEN t.strict_lag2_confirmation_date IS NOT NULL THEN 1 ELSE 0 END
      AS meningioma_event_confirmed,

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

FROM revised_eligible_base a
JOIN meningioma_temporal_summary_revised t USING (patient_id)
WHERE t.broad_preindex_date IS NULL
  AND t.broad_lag2_date IS NULL
  AND a.last_encounter_date > a.index_date + INTERVAL 730 DAY;
")

# ------------------------------------------------------------
# 6. QC
# ------------------------------------------------------------
revised_counts <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    COUNT(*) AS patients,
    SUM(meningioma_event_strict) AS strict_events,
    SUM(meningioma_event_confirmed) AS confirmed_events,
    SUM(followup_days_strict)/365.25 AS total_person_years,
    100000.0*SUM(meningioma_event_strict)/(SUM(followup_days_strict)/365.25)
      AS crude_incidence_per_100k_py,
    MEDIAN(followup_days_strict)/365.25 AS median_followup_y
FROM analysis_cohort_primary_revised
GROUP BY exposure
ORDER BY exposure;
"))

lag2_counts <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    COUNT(*) AS patients,
    SUM(meningioma_event_strict) AS strict_events,
    SUM(meningioma_event_confirmed) AS confirmed_events,
    SUM(followup_days_strict)/365.25 AS total_person_years,
    100000.0*SUM(meningioma_event_strict)/(SUM(followup_days_strict)/365.25)
      AS crude_incidence_per_100k_py,
    MEDIAN(followup_days_strict)/365.25 AS median_followup_y
FROM analysis_cohort_lag2_revised
GROUP BY exposure
ORDER BY exposure;
"))

compare_old_new <- as.data.table(dbGetQuery(con, "
WITH old AS (
  SELECT exposure, COUNT(*) AS old_n
  FROM analysis_cohort_primary_corrected
  GROUP BY exposure
),
new AS (
  SELECT exposure, COUNT(*) AS revised_n
  FROM analysis_cohort_primary_revised
  GROUP BY exposure
)
SELECT
  new.exposure,
  old.old_n,
  new.revised_n,
  new.revised_n - old.old_n AS patients_added
FROM new
LEFT JOIN old USING (exposure)
ORDER BY exposure;
"))

fwrite(revised_counts, file.path(QC_DIR, "05_revised_primary_counts.csv"))
fwrite(lag2_counts, file.path(QC_DIR, "05_revised_lag2_counts.csv"))
fwrite(compare_old_new, file.path(QC_DIR, "05_old_vs_revised_counts.csv"))

# ------------------------------------------------------------
# 7. EXPORT
# ------------------------------------------------------------
p1 <- sql_path(file.path(PARQUET_DIR, "analysis_cohort_primary_revised.parquet"))
p2 <- sql_path(file.path(PARQUET_DIR, "analysis_cohort_lag2_revised.parquet"))

dbExecute(con, sprintf(
  "COPY analysis_cohort_primary_revised TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000);", p1
))
dbExecute(con, sprintf(
  "COPY analysis_cohort_lag2_revised TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000);", p2
))

cat("\n================ OLD VS REVISED ====================\n")
print(compare_old_new)

cat("\n================ REVISED PRIMARY ===================\n")
print(revised_counts)

cat("\n================ REVISED 2-YEAR LAG ================\n")
print(lag2_counts)

cat("\nStage 5 cohort rebuild complete.\n")

dbDisconnect(con, shutdown = TRUE)
