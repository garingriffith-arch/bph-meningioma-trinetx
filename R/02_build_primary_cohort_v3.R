# ============================================================
# 02_build_primary_cohort_v3.R
# Finasteride / Dutasteride and Incident Meningioma Study
#
# PATCHED VERSION: avoids base-R sprintf() 8192-byte limit and can resume
# after a failed prior run without repeating the completed encounter scan.
#
# Builds the primary active-comparator, new-user cohort:
#   5-alpha reductase inhibitor (5ARI) initiation
#       vs
#   uroselective alpha-blocker initiation
#
# Requires:
#   00_validate_export.R
#   01_extract_medication_timeline.R
#
# Heavy scans in this script:
#   encounter.csv  : one pass
#   diagnosis.csv  : one pass
#
# Raw CSVs are NEVER modified.
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

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------
RAW_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/Data/BPH"
R_DIR   <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"

WORK_DIR <- "C:/TriNetX/bph_m_work"
DB_DIR <- file.path(WORK_DIR, "duckdb")
TEMP_DIR <- file.path(WORK_DIR, "duckdb_temp")
QC_DIR <- file.path(WORK_DIR, "qc")
PARQUET_DIR <- file.path(WORK_DIR, "parquet")

dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PARQUET_DIR, recursive = TRUE, showWarnings = FALSE)

DB_PATH <- file.path(DB_DIR, "bph_m.duckdb")

ENCOUNTER_FILE <- file.path(RAW_DIR, "encounter.csv")
DIAGNOSIS_FILE <- file.path(RAW_DIR, "diagnosis.csv")

if (!file.exists(ENCOUNTER_FILE)) stop("Missing: ", ENCOUNTER_FILE)
if (!file.exists(DIAGNOSIS_FILE)) stop("Missing: ", DIAGNOSIS_FILE)

# ------------------------------------------------------------
# 2. STUDY PARAMETERS
# ------------------------------------------------------------
# Calendar floor:
# Removes implausible legacy/default dates and focuses on the modern EHR era.
INDEX_START <- as.Date("2004-01-01")

# Baseline requirements
MIN_AGE <- 40L
MAX_AGE <- 89L

LOOKBACK_DAYS <- 365L
BASELINE_COVARIATE_DAYS <- 1095L    # 3 years

# Require a BPH diagnosis around treatment initiation.
BPH_LOOKBACK_DAYS <- 365L
BPH_GRACE_DAYS <- 30L

# Exposure stabilization period:
# Exclude patients who add the opposite study class within 30 days.
TREATMENT_GRACE_DAYS <- 30L

# Primary outcome lag:
# meningiomas during the first year are not counted as incident outcomes.
OUTCOME_LAG_DAYS <- 365L

# We require at least 365 days of possible outcome observation AFTER the lag.
# Therefore index must be at least 730 days before the dataset's observed cutoff.
MIN_POST_LAG_DAYS <- 365L

# Known finasteride 1-mg initiators are excluded from the PRIMARY cohort because
# 1 mg is strongly associated with the alopecia indication.
EXCLUDE_KNOWN_FINASTERIDE_1MG <- TRUE

# If a previous run already completed an expensive raw-file scan before
# stopping, reuse the materialized DuckDB table rather than scanning again.
# This is especially useful after the original sprintf() failure, because
# candidate_observation was usually completed before the error occurred.
REUSE_EXISTING_HEAVY_TABLES <- TRUE

# ------------------------------------------------------------
# 3. PERFORMANCE
# ------------------------------------------------------------
N_THREADS <- min(24L, max(4L, parallel::detectCores(logical = TRUE) - 6L))
MEMORY_LIMIT <- "48GB"

sql_path <- function(x) {
  x <- normalizePath(x, winslash = "/", mustWork = FALSE)
  gsub("'", "''", x, fixed = TRUE)
}

sql_date <- function(x) format(as.Date(x), "%Y-%m-%d")

con <- dbConnect(duckdb::duckdb(), dbdir = DB_PATH, read_only = FALSE)

on.exit({
  try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

dbExecute(con, sprintf("PRAGMA threads=%d;", N_THREADS))
dbExecute(con, sprintf("PRAGMA memory_limit='%s';", MEMORY_LIMIT))
dbExecute(con, sprintf("PRAGMA temp_directory='%s';", sql_path(TEMP_DIR)))
dbExecute(con, "SET preserve_insertion_order=false;")

needed <- c(
  "patient_core",
  "study_med_events",
  "study_med_firsts",
  "index_candidates_prelim",
  "class_sequence"
)
missing <- setdiff(needed, dbListTables(con))
if (length(missing) > 0) {
  stop(
    "Missing prior-stage DuckDB tables: ",
    paste(missing, collapse = ", "),
    "\nRun Stage 1 first."
  )
}

# ------------------------------------------------------------
# 4. BUILD CLEAN PRE-CLINICAL-HISTORY INDEX CANDIDATES
# ------------------------------------------------------------
# Primary comparator:
#   5ARI vs uroselective alpha blocker
#
# We require one study drug on the first study-medication day so that the
# initiating strategy is unambiguous.
#
# Age is calculated from year of birth because the export provides YOB rather
# than an exact DOB.
dbExecute(con, "
CREATE OR REPLACE TABLE candidate_index0 AS
SELECT
    i.patient_id,
    i.index_date_prelim AS index_date,
    i.exposure_prelim AS exposure,
    i.first_day_drugs AS index_drug,
    i.n_drugs_on_first_day,

    p.sex,
    p.race,
    p.ethnicity,
    p.marital_status,
    p.year_of_birth,
    p.patient_regional_location,
    p.source_id AS patient_source_id,

    YEAR(i.index_date_prelim) - p.year_of_birth AS age_at_index,

    i.first_finasteride,
    i.first_dutasteride,
    i.first_tamsulosin,
    i.first_alfuzosin,
    i.first_silodosin,
    i.first_doxazosin,
    i.first_terazosin,
    i.first_5ari,
    i.first_uro_alpha,
    i.first_nonuro_alpha

FROM index_candidates_prelim i
JOIN patient_core p USING (patient_id)

WHERE i.exposure_prelim IN ('5ARI', 'URO_ALPHA')
  AND i.n_drugs_on_first_day = 1;
")

try(dbExecute(con,
  "CREATE INDEX IF NOT EXISTS idx_candidate_index0_patient
   ON candidate_index0(patient_id);"
), silent = TRUE)

# ------------------------------------------------------------
# 5. IDENTIFY INDEX-DAY 5ARI STRENGTH
# ------------------------------------------------------------
# This lets us remove known finasteride 1-mg initiators from the primary
# analysis while retaining patients with missing strength for sensitivity work.
dbExecute(con, "
CREATE OR REPLACE TABLE candidate_index_strength AS
SELECT
    c.patient_id,

    MIN(
      CASE
        WHEN c.exposure = '5ARI'
         AND e.start_date = c.index_date
         AND e.drug = 'finasteride'
        THEN e.strength
      END
    ) AS index_finasteride_strength,

    MIN(
      CASE
        WHEN c.exposure = '5ARI'
         AND e.start_date = c.index_date
         AND e.drug = 'dutasteride'
        THEN e.strength
      END
    ) AS index_dutasteride_strength

FROM candidate_index0 c
LEFT JOIN study_med_events e
  ON e.patient_id = c.patient_id
 AND e.start_date = c.index_date

GROUP BY c.patient_id;
")

# ------------------------------------------------------------
# 6. SCAN ENCOUNTER.CSV ONCE
# ------------------------------------------------------------
# We obtain:
#   - first/last observable encounter
#   - activity during the 365-day lookback
#   - utilization during 1-year and 3-year baseline periods
#   - ED / inpatient / ambulatory utilization
#
# This does NOT claim continuous insurance enrollment.
if (REUSE_EXISTING_HEAVY_TABLES &&
    "candidate_observation" %in% dbListTables(con)) {

  cat("\n============================================================\n")
  cat("Reusing existing candidate_observation table; encounter.csv will NOT be rescanned.\n")
  cat("============================================================\n")

} else {

  cat("\n============================================================\n")
  cat("Scanning encounter.csv once...\n")
  cat("============================================================\n")
  enc_t0 <- Sys.time()

  dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE candidate_observation AS
  WITH enc AS (
      SELECT
          patient_id,
          TRY_CAST(STRPTIME(start_date, '%%Y%%m%%d') AS DATE) AS encounter_date,
          type AS encounter_type
      FROM read_csv_auto(
          '%s',
          header=true,
          all_varchar=true,
          ignore_errors=false
      )
      WHERE start_date IS NOT NULL
        AND start_date <> ''
  )
  SELECT
      c.patient_id,

      MIN(e.encounter_date) AS first_encounter_date,
      MAX(e.encounter_date) AS last_encounter_date,

      COUNT(*) FILTER (
          WHERE e.encounter_date >= c.index_date - INTERVAL %d DAY
            AND e.encounter_date <  c.index_date
      ) AS n_encounters_prior_365,

      COUNT(*) FILTER (
          WHERE e.encounter_date >= c.index_date - INTERVAL %d DAY
            AND e.encounter_date <  c.index_date
      ) AS n_encounters_prior_1095,

      COUNT(*) FILTER (
          WHERE e.encounter_date >= c.index_date - INTERVAL %d DAY
            AND e.encounter_date <  c.index_date
            AND e.encounter_type = 'AMB'
      ) AS n_ambulatory_prior_365,

      COUNT(*) FILTER (
          WHERE e.encounter_date >= c.index_date - INTERVAL %d DAY
            AND e.encounter_date <  c.index_date
            AND e.encounter_type IN ('EMER')
      ) AS n_ed_prior_365,

      COUNT(*) FILTER (
          WHERE e.encounter_date >= c.index_date - INTERVAL %d DAY
            AND e.encounter_date <  c.index_date
            AND e.encounter_type IN ('IMP','ACUTE','NONAC','SS')
      ) AS n_inpatient_prior_365

  FROM candidate_index0 c
  LEFT JOIN enc e
    ON e.patient_id = c.patient_id

  GROUP BY c.patient_id;
  ",
  sql_path(ENCOUNTER_FILE),
  LOOKBACK_DAYS,
  BASELINE_COVARIATE_DAYS,
  LOOKBACK_DAYS,
  LOOKBACK_DAYS,
  LOOKBACK_DAYS
  ))

  cat("encounter.csv scan completed in ",
      round(as.numeric(difftime(Sys.time(), enc_t0, units = "mins")), 1),
      " minutes.\n", sep = "")

}

# Dataset observational cutoff based on actual encounters.
data_cutoff <- dbGetQuery(
  con,
  "SELECT MAX(last_encounter_date) AS data_cutoff FROM candidate_observation;"
)$data_cutoff[1]

data_cutoff <- as.Date(data_cutoff)

if (is.na(data_cutoff)) stop("Could not determine data cutoff from encounter.csv")

INDEX_END <- data_cutoff - (OUTCOME_LAG_DAYS + MIN_POST_LAG_DAYS)

cat("\nObserved data cutoff:", format(data_cutoff), "\n")
cat("Primary index window:", format(INDEX_START), "through",
    format(INDEX_END), "\n")

# ------------------------------------------------------------
# 7. SCAN DIAGNOSIS.CSV ONCE
# ------------------------------------------------------------
# Code normalization:
#   - uppercase
#   - remove decimal point
#
# Phenotypes retained:
#
# Required cohort definitions:
#   BPH
#   prior prostate cancer
#   prior meningioma
#
# Outcomes:
#   STRICT intracranial meningioma:
#       ICD-10-CM D32.0, D32.9
#       ICD-9-CM 225.2
#
#   BROAD meningeal tumor sensitivity:
#       strict +
#       D42.0 / D42.9 (uncertain behavior)
#       C70.0 / C70.9 (malignant meninges)
#       ICD-9 237.6 / 192.1
#
# Baseline confounders (3-y pre-index):
#   hypertension, diabetes, obesity, CKD, nicotine/smoking, CAD, HF,
#   cerebrovascular disease, hyperlipidemia, osteoporosis,
#   androgenic alopecia, urinary retention, UTI, hematuria, bladder stone
if (REUSE_EXISTING_HEAVY_TABLES &&
    "candidate_dx_summary" %in% dbListTables(con)) {

  cat("\n============================================================\n")
  cat("Reusing existing candidate_dx_summary table; diagnosis.csv will NOT be rescanned.\n")
  cat("============================================================\n")

} else {

  cat("\n============================================================\n")
  cat("Scanning diagnosis.csv once...\n")
  cat("============================================================\n")
  dx_t0 <- Sys.time()

  # IMPORTANT:
  # Do NOT put this long SQL statement inside sprintf().
  # Base R sprintf() has an 8192-byte maximum format-string length.
  # The diagnosis query is longer than that, which caused:
  #   'fmt' length exceeds maximal format length 8192
  #
  # Instead we use explicit {{TOKEN}} placeholders and fixed-string
  # substitution. The SQL itself can be much longer than 8192 bytes.

  sql_fill <- function(template, values) {
    out <- template
    for (nm in names(values)) {
      token <- paste0("{{", nm, "}}")
      out <- gsub(
        token,
        as.character(values[[nm]]),
        out,
        fixed = TRUE
      )
    }
    out
  }

  dx_sql_template <- "
CREATE OR REPLACE TABLE candidate_dx_summary AS

WITH dxraw AS (
    SELECT
        patient_id,
        UPPER(TRIM(code_system)) AS code_system,
        UPPER(REPLACE(TRIM(code), '.', '')) AS code_norm,
        TRY_CAST(STRPTIME(date, '%Y%m%d') AS DATE) AS dx_date
    FROM read_csv_auto(
        '{{DIAGNOSIS_FILE}}',
        header=true,
        all_varchar=true,
        ignore_errors=false
    )
    WHERE date IS NOT NULL
      AND date <> ''
      AND code IS NOT NULL
      AND code <> ''
      AND UPPER(TRIM(code_system)) IN ('ICD-10-CM', 'ICD-9-CM')
),

-- Keep only diagnosis families actually used by this study before the
-- patient-index join. This still scans diagnosis.csv once, but it avoids
-- carrying tens/hundreds of millions of irrelevant diagnoses downstream.
dx0 AS (
    SELECT *
    FROM dxraw
    WHERE
      (
        code_system = 'ICD-10-CM'
        AND (
             code_norm LIKE 'N40%'
          OR code_norm = 'C61'
          OR code_norm IN ('D320','D329','D420','D429','C700','C709')
          OR code_norm LIKE 'I10%'
          OR code_norm LIKE 'I11%'
          OR code_norm LIKE 'I12%'
          OR code_norm LIKE 'I13%'
          OR code_norm LIKE 'I15%'
          OR code_norm LIKE 'I16%'
          OR code_norm LIKE 'E08%'
          OR code_norm LIKE 'E09%'
          OR code_norm LIKE 'E10%'
          OR code_norm LIKE 'E11%'
          OR code_norm LIKE 'E12%'
          OR code_norm LIKE 'E13%'
          OR code_norm LIKE 'E66%'
          OR code_norm LIKE 'N18%'
          OR code_norm LIKE 'F17%'
          OR code_norm IN ('Z720','Z87891')
          OR code_norm LIKE 'I25%'
          OR code_norm LIKE 'I50%'
          OR (code_norm >= 'I60' AND code_norm < 'I70')
          OR code_norm LIKE 'E78%'
          OR code_norm LIKE 'M80%'
          OR code_norm LIKE 'M81%'
          OR code_norm LIKE 'L64%'
          OR code_norm LIKE 'R33%'
          OR code_norm = 'N390'
          OR code_norm LIKE 'R31%'
          OR code_norm LIKE 'N21%'
        )
      )
      OR
      (
        code_system = 'ICD-9-CM'
        AND (
             code_norm LIKE '600%'
          OR code_norm = '185'
          OR code_norm IN ('2252','2376','1921')
          OR code_norm LIKE '401%'
          OR code_norm LIKE '402%'
          OR code_norm LIKE '403%'
          OR code_norm LIKE '404%'
          OR code_norm LIKE '405%'
          OR code_norm LIKE '250%'
          OR code_norm LIKE '2780%'
          OR code_norm = '278'
          OR code_norm LIKE '585%'
          OR code_norm LIKE '3051%'
          OR code_norm = 'V1582'
          OR code_norm LIKE '414%'
          OR code_norm LIKE '428%'
          OR (code_norm >= '430' AND code_norm < '439')
          OR code_norm LIKE '272%'
          OR code_norm LIKE '7330%'
          OR code_norm LIKE '704%'
          OR code_norm LIKE '7882%'
          OR code_norm = '5990'
          OR code_norm LIKE '5997%'
          OR code_norm LIKE '594%'
        )
      )
),

dx AS (
    SELECT
        d.*,

        CASE
          WHEN code_system = 'ICD-10-CM'
               AND code_norm LIKE 'N40%'
            THEN 1
          WHEN code_system = 'ICD-9-CM'
               AND code_norm LIKE '600%'
            THEN 1
          ELSE 0
        END AS is_bph,

        CASE
          WHEN code_system = 'ICD-10-CM'
               AND code_norm = 'C61'
            THEN 1
          WHEN code_system = 'ICD-9-CM'
               AND code_norm = '185'
            THEN 1
          ELSE 0
        END AS is_prostate_cancer,

        CASE
          WHEN code_system = 'ICD-10-CM'
               AND code_norm IN ('D320','D329')
            THEN 1
          WHEN code_system = 'ICD-9-CM'
               AND code_norm = '2252'
            THEN 1
          ELSE 0
        END AS is_meningioma_strict,

        CASE
          WHEN code_system = 'ICD-10-CM'
               AND code_norm IN ('D320','D329','D420','D429','C700','C709')
            THEN 1
          WHEN code_system = 'ICD-9-CM'
               AND code_norm IN ('2252','2376','1921')
            THEN 1
          ELSE 0
        END AS is_meningeal_tumor_broad,

        CASE
          WHEN (code_system='ICD-10-CM' AND
                (code_norm LIKE 'I10%' OR code_norm LIKE 'I11%' OR
                 code_norm LIKE 'I12%' OR code_norm LIKE 'I13%' OR
                 code_norm LIKE 'I15%' OR code_norm LIKE 'I16%'))
            OR (code_system='ICD-9-CM' AND
                (code_norm LIKE '401%' OR code_norm LIKE '402%' OR
                 code_norm LIKE '403%' OR code_norm LIKE '404%' OR
                 code_norm LIKE '405%'))
          THEN 1 ELSE 0
        END AS is_hypertension,

        CASE
          WHEN (code_system='ICD-10-CM' AND
                (code_norm LIKE 'E08%' OR code_norm LIKE 'E09%' OR
                 code_norm LIKE 'E10%' OR code_norm LIKE 'E11%' OR
                 code_norm LIKE 'E12%' OR code_norm LIKE 'E13%'))
            OR (code_system='ICD-9-CM' AND code_norm LIKE '250%')
          THEN 1 ELSE 0
        END AS is_diabetes,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'E66%')
            OR (code_system='ICD-9-CM' AND
                (code_norm LIKE '2780%' OR code_norm='278'))
          THEN 1 ELSE 0
        END AS is_obesity,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'N18%')
            OR (code_system='ICD-9-CM' AND code_norm LIKE '585%')
          THEN 1 ELSE 0
        END AS is_ckd,

        CASE
          WHEN (code_system='ICD-10-CM' AND
                (code_norm LIKE 'F17%' OR code_norm='Z720' OR
                 code_norm='Z87891'))
            OR (code_system='ICD-9-CM' AND
                (code_norm LIKE '3051%' OR code_norm='V1582'))
          THEN 1 ELSE 0
        END AS is_smoking,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'I25%')
            OR (code_system='ICD-9-CM' AND code_norm LIKE '414%')
          THEN 1 ELSE 0
        END AS is_cad,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'I50%')
            OR (code_system='ICD-9-CM' AND code_norm LIKE '428%')
          THEN 1 ELSE 0
        END AS is_hf,

        CASE
          WHEN (code_system='ICD-10-CM' AND
                code_norm >= 'I60' AND code_norm < 'I70')
            OR (code_system='ICD-9-CM' AND
                code_norm >= '430' AND code_norm < '439')
          THEN 1 ELSE 0
        END AS is_cerebrovascular,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'E78%')
            OR (code_system='ICD-9-CM' AND code_norm LIKE '272%')
          THEN 1 ELSE 0
        END AS is_hyperlipidemia,

        CASE
          WHEN (code_system='ICD-10-CM' AND
                (code_norm LIKE 'M80%' OR code_norm LIKE 'M81%'))
            OR (code_system='ICD-9-CM' AND code_norm LIKE '7330%')
          THEN 1 ELSE 0
        END AS is_osteoporosis,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'L64%')
            OR (code_system='ICD-9-CM' AND code_norm LIKE '704%')
          THEN 1 ELSE 0
        END AS is_alopecia,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'R33%')
            OR (code_system='ICD-9-CM' AND code_norm LIKE '7882%')
          THEN 1 ELSE 0
        END AS is_urinary_retention,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm='N390')
            OR (code_system='ICD-9-CM' AND code_norm='5990')
          THEN 1 ELSE 0
        END AS is_uti,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'R31%')
            OR (code_system='ICD-9-CM' AND code_norm LIKE '5997%')
          THEN 1 ELSE 0
        END AS is_hematuria,

        CASE
          WHEN (code_system='ICD-10-CM' AND code_norm LIKE 'N21%')
            OR (code_system='ICD-9-CM' AND code_norm LIKE '594%')
          THEN 1 ELSE 0
        END AS is_bladder_stone

    FROM dx0 d
),

j AS (
    SELECT
        c.patient_id,
        c.index_date,
        d.*
    FROM candidate_index0 c
    LEFT JOIN dx d
      ON d.patient_id = c.patient_id
)

SELECT
    patient_id,

    MIN(dx_date) FILTER (WHERE is_bph = 1) AS first_bph_date,
    MAX(dx_date) FILTER (WHERE is_bph = 1) AS last_bph_date,

    MAX(CASE
        WHEN is_bph = 1
         AND dx_date >= index_date - INTERVAL {{BPH_LOOKBACK_DAYS}} DAY
         AND dx_date <= index_date + INTERVAL {{BPH_GRACE_DAYS}} DAY
        THEN 1 ELSE 0 END
    ) AS bph_near_index,

    MAX(CASE
        WHEN is_prostate_cancer = 1
         AND dx_date < index_date
        THEN 1 ELSE 0 END
    ) AS prior_prostate_cancer,

    MAX(CASE
        WHEN is_meningioma_strict = 1
         AND dx_date < index_date
        THEN 1 ELSE 0 END
    ) AS prior_meningioma_strict,

    MAX(CASE
        WHEN is_meningeal_tumor_broad = 1
         AND dx_date < index_date
        THEN 1 ELSE 0 END
    ) AS prior_meningeal_tumor_broad,

    MIN(dx_date) FILTER (
        WHERE is_meningioma_strict = 1
          AND dx_date >= index_date + INTERVAL {{OUTCOME_LAG_DAYS}} DAY
    ) AS incident_meningioma_date_strict,

    MIN(dx_date) FILTER (
        WHERE is_meningeal_tumor_broad = 1
          AND dx_date >= index_date + INTERVAL {{OUTCOME_LAG_DAYS}} DAY
    ) AS incident_meningeal_tumor_date_broad,

    MAX(CASE WHEN is_hypertension=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS hypertension,

    MAX(CASE WHEN is_diabetes=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS diabetes,

    MAX(CASE WHEN is_obesity=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS obesity,

    MAX(CASE WHEN is_ckd=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS ckd,

    MAX(CASE WHEN is_smoking=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS smoking,

    MAX(CASE WHEN is_cad=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS cad,

    MAX(CASE WHEN is_hf=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS heart_failure,

    MAX(CASE WHEN is_cerebrovascular=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS cerebrovascular_disease,

    MAX(CASE WHEN is_hyperlipidemia=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS hyperlipidemia,

    MAX(CASE WHEN is_osteoporosis=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS osteoporosis,

    MAX(CASE WHEN is_alopecia=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS alopecia,

    MAX(CASE WHEN is_urinary_retention=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS urinary_retention,

    MAX(CASE WHEN is_uti=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS uti,

    MAX(CASE WHEN is_hematuria=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS hematuria,

    MAX(CASE WHEN is_bladder_stone=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_COVARIATE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS bladder_stone

FROM j
GROUP BY patient_id;
"

  dx_sql <- sql_fill(
    dx_sql_template,
    list(
      DIAGNOSIS_FILE = sql_path(DIAGNOSIS_FILE),
      BPH_LOOKBACK_DAYS = BPH_LOOKBACK_DAYS,
      BPH_GRACE_DAYS = BPH_GRACE_DAYS,
      OUTCOME_LAG_DAYS = OUTCOME_LAG_DAYS,
      BASELINE_COVARIATE_DAYS = BASELINE_COVARIATE_DAYS
    )
  )

  dbExecute(con, dx_sql)

  cat("diagnosis.csv scan completed in ",
      round(as.numeric(difftime(Sys.time(), dx_t0, units = "mins")), 1),
      " minutes.\n", sep = "")

}

# ------------------------------------------------------------
# 8. MERGE EVERYTHING + ELIGIBILITY FLAGS
# ------------------------------------------------------------
dbExecute(con, sprintf("
CREATE OR REPLACE TABLE cohort_eligibility AS
SELECT
    c.*,

    s.index_finasteride_strength,
    s.index_dutasteride_strength,

    o.first_encounter_date,
    o.last_encounter_date,
    o.n_encounters_prior_365,
    o.n_encounters_prior_1095,
    o.n_ambulatory_prior_365,
    o.n_ed_prior_365,
    o.n_inpatient_prior_365,

    d.first_bph_date,
    d.last_bph_date,
    COALESCE(d.bph_near_index, 0) AS bph_near_index,
    COALESCE(d.prior_prostate_cancer, 0) AS prior_prostate_cancer,
    COALESCE(d.prior_meningioma_strict, 0) AS prior_meningioma_strict,
    COALESCE(d.prior_meningeal_tumor_broad, 0) AS prior_meningeal_tumor_broad,

    d.incident_meningioma_date_strict,
    d.incident_meningeal_tumor_date_broad,

    COALESCE(d.hypertension,0) AS hypertension,
    COALESCE(d.diabetes,0) AS diabetes,
    COALESCE(d.obesity,0) AS obesity,
    COALESCE(d.ckd,0) AS ckd,
    COALESCE(d.smoking,0) AS smoking,
    COALESCE(d.cad,0) AS cad,
    COALESCE(d.heart_failure,0) AS heart_failure,
    COALESCE(d.cerebrovascular_disease,0) AS cerebrovascular_disease,
    COALESCE(d.hyperlipidemia,0) AS hyperlipidemia,
    COALESCE(d.osteoporosis,0) AS osteoporosis,
    COALESCE(d.alopecia,0) AS alopecia,
    COALESCE(d.urinary_retention,0) AS urinary_retention,
    COALESCE(d.uti,0) AS uti,
    COALESCE(d.hematuria,0) AS hematuria,
    COALESCE(d.bladder_stone,0) AS bladder_stone,

    CASE
      WHEN c.exposure = '5ARI'
       AND c.first_uro_alpha IS NOT NULL
       AND c.first_uro_alpha > c.index_date
       AND c.first_uro_alpha <= c.index_date + INTERVAL %d DAY
      THEN 1

      WHEN c.exposure = 'URO_ALPHA'
       AND c.first_5ari IS NOT NULL
       AND c.first_5ari > c.index_date
       AND c.first_5ari <= c.index_date + INTERVAL %d DAY
      THEN 1

      ELSE 0
    END AS opposite_class_within_30d,

    CASE
      WHEN c.sex = 'M' THEN 1 ELSE 0
    END AS eligible_male,

    CASE
      WHEN c.age_at_index BETWEEN %d AND %d THEN 1 ELSE 0
    END AS eligible_age,

    CASE
      WHEN c.index_date >= DATE '%s'
       AND c.index_date <= DATE '%s'
      THEN 1 ELSE 0
    END AS eligible_calendar,

    CASE
      WHEN o.first_encounter_date <= c.index_date - INTERVAL %d DAY
       AND o.n_encounters_prior_365 >= 1
      THEN 1 ELSE 0
    END AS eligible_lookback,

    CASE
      WHEN o.last_encounter_date >= c.index_date + INTERVAL %d DAY
      THEN 1 ELSE 0
    END AS eligible_postlag_observation,

    CASE
      WHEN COALESCE(d.bph_near_index,0) = 1 THEN 1 ELSE 0
    END AS eligible_bph,

    CASE
      WHEN COALESCE(d.prior_prostate_cancer,0) = 0 THEN 1 ELSE 0
    END AS eligible_no_prostate_cancer,

    CASE
      WHEN COALESCE(d.prior_meningeal_tumor_broad,0) = 0 THEN 1 ELSE 0
    END AS eligible_no_prior_meningeal_tumor,

    CASE
      WHEN NOT (
            c.exposure = '5ARI'
        AND c.index_drug = 'finasteride'
        AND UPPER(TRIM(COALESCE(s.index_finasteride_strength,''))) = '1 MG'
      )
      THEN 1 ELSE 0
    END AS eligible_not_known_finasteride_1mg

FROM candidate_index0 c
LEFT JOIN candidate_index_strength s USING (patient_id)
LEFT JOIN candidate_observation o USING (patient_id)
LEFT JOIN candidate_dx_summary d USING (patient_id);
",
TREATMENT_GRACE_DAYS,
TREATMENT_GRACE_DAYS,
MIN_AGE,
MAX_AGE,
sql_date(INDEX_START),
sql_date(INDEX_END),
LOOKBACK_DAYS,
OUTCOME_LAG_DAYS + MIN_POST_LAG_DAYS
))

# ------------------------------------------------------------
# 9. COHORT FLOW
# ------------------------------------------------------------
# Build the named conditions programmatically.
# In R, an expression such as:
#   sprintf("Age %d-%d", MIN_AGE, MAX_AGE) =
# is invalid syntax on the left side of "=" inside c()/list().
#
# Using named assignment avoids that parser error and is easier to maintain.

flow_queries <- character(0)

flow_queries["Initial unambiguous 5ARI/URO-alpha initiators"] <-
  "TRUE"

flow_queries["Male"] <-
  "eligible_male = 1"

flow_queries[sprintf("Age %d-%d", MIN_AGE, MAX_AGE)] <-
  "eligible_male = 1 AND eligible_age = 1"

flow_queries[sprintf("Index %s to %s", INDEX_START, INDEX_END)] <-
  paste0(
    "eligible_male = 1 ",
    "AND eligible_age = 1 ",
    "AND eligible_calendar = 1"
  )

flow_queries[">=365 d prior observable history + encounter in prior year"] <-
  paste0(
    "eligible_male = 1 ",
    "AND eligible_age = 1 ",
    "AND eligible_calendar = 1 ",
    "AND eligible_lookback = 1"
  )

flow_queries[
  sprintf(
    "BPH diagnosis within %d d before to %d d after index",
    BPH_LOOKBACK_DAYS,
    BPH_GRACE_DAYS
  )
] <-
  paste0(
    "eligible_male = 1 ",
    "AND eligible_age = 1 ",
    "AND eligible_calendar = 1 ",
    "AND eligible_lookback = 1 ",
    "AND eligible_bph = 1"
  )

flow_queries["No prior prostate cancer"] <-
  paste0(
    "eligible_male = 1 ",
    "AND eligible_age = 1 ",
    "AND eligible_calendar = 1 ",
    "AND eligible_lookback = 1 ",
    "AND eligible_bph = 1 ",
    "AND eligible_no_prostate_cancer = 1"
  )

flow_queries["No prior meningeal tumor"] <-
  paste0(
    "eligible_male = 1 ",
    "AND eligible_age = 1 ",
    "AND eligible_calendar = 1 ",
    "AND eligible_lookback = 1 ",
    "AND eligible_bph = 1 ",
    "AND eligible_no_prostate_cancer = 1 ",
    "AND eligible_no_prior_meningeal_tumor = 1"
  )

flow_queries["No opposite study class added within 30 d"] <-
  paste0(
    "eligible_male = 1 ",
    "AND eligible_age = 1 ",
    "AND eligible_calendar = 1 ",
    "AND eligible_lookback = 1 ",
    "AND eligible_bph = 1 ",
    "AND eligible_no_prostate_cancer = 1 ",
    "AND eligible_no_prior_meningeal_tumor = 1 ",
    "AND opposite_class_within_30d = 0"
  )

flow_queries["Exclude known finasteride 1-mg initiators"] <-
  paste0(
    "eligible_male = 1 ",
    "AND eligible_age = 1 ",
    "AND eligible_calendar = 1 ",
    "AND eligible_lookback = 1 ",
    "AND eligible_bph = 1 ",
    "AND eligible_no_prostate_cancer = 1 ",
    "AND eligible_no_prior_meningeal_tumor = 1 ",
    "AND opposite_class_within_30d = 0 ",
    if (EXCLUDE_KNOWN_FINASTERIDE_1MG)
      "AND eligible_not_known_finasteride_1mg = 1"
    else
      ""
  )

flow_queries[
  sprintf(
    ">=%d d observable follow-up after %d-d outcome lag",
    MIN_POST_LAG_DAYS,
    OUTCOME_LAG_DAYS
  )
] <-
  paste0(
    "eligible_male = 1 ",
    "AND eligible_age = 1 ",
    "AND eligible_calendar = 1 ",
    "AND eligible_lookback = 1 ",
    "AND eligible_bph = 1 ",
    "AND eligible_no_prostate_cancer = 1 ",
    "AND eligible_no_prior_meningeal_tumor = 1 ",
    "AND opposite_class_within_30d = 0 ",
    if (EXCLUDE_KNOWN_FINASTERIDE_1MG)
      "AND eligible_not_known_finasteride_1mg = 1 "
    else
      "",
    "AND eligible_postlag_observation = 1"
  )

flow <- rbindlist(lapply(seq_along(flow_queries), function(k) {
  criterion <- names(flow_queries)[k]
  condition <- flow_queries[[k]]

  x <- dbGetQuery(
    con,
    sprintf("
      SELECT
        COUNT(*) AS total,
        SUM(CASE WHEN exposure='5ARI' THEN 1 ELSE 0 END) AS n_5ari,
        SUM(CASE WHEN exposure='URO_ALPHA' THEN 1 ELSE 0 END) AS n_uro_alpha
      FROM cohort_eligibility
      WHERE %s;
    ", condition)
  )

  data.table(
    step = k,
    criterion = criterion,
    total = as.numeric(x$total),
    n_5ari = as.numeric(x$n_5ari),
    n_uro_alpha = as.numeric(x$n_uro_alpha)
  )
}))

flow[, excluded_since_prior := shift(total, 1L, type = "lag") - total]
flow[1, excluded_since_prior := NA_real_]

fwrite(flow, file.path(QC_DIR, "02_cohort_flow.csv"))

# ------------------------------------------------------------
# 10. CREATE PRIMARY COHORT
# ------------------------------------------------------------
primary_where <- paste0(
  "eligible_male = 1 ",
  "AND eligible_age = 1 ",
  "AND eligible_calendar = 1 ",
  "AND eligible_lookback = 1 ",
  "AND eligible_bph = 1 ",
  "AND eligible_no_prostate_cancer = 1 ",
  "AND eligible_no_prior_meningeal_tumor = 1 ",
  "AND opposite_class_within_30d = 0 ",
  if (EXCLUDE_KNOWN_FINASTERIDE_1MG)
    "AND eligible_not_known_finasteride_1mg = 1 "
  else "",
  "AND eligible_postlag_observation = 1"
)

dbExecute(con, sprintf("
CREATE OR REPLACE TABLE analysis_cohort_primary AS
SELECT
    *,

    index_date + INTERVAL %d DAY AS landmark_date,

    CASE
      WHEN incident_meningioma_date_strict IS NOT NULL
       AND incident_meningioma_date_strict <= last_encounter_date
      THEN 1 ELSE 0
    END AS meningioma_event_strict,

    CASE
      WHEN incident_meningeal_tumor_date_broad IS NOT NULL
       AND incident_meningeal_tumor_date_broad <= last_encounter_date
      THEN 1 ELSE 0
    END AS meningeal_event_broad,

    CASE
      WHEN incident_meningioma_date_strict IS NOT NULL
       AND incident_meningioma_date_strict <= last_encounter_date
      THEN incident_meningioma_date_strict
      ELSE last_encounter_date
    END AS followup_end_strict,

    CASE
      WHEN incident_meningeal_tumor_date_broad IS NOT NULL
       AND incident_meningeal_tumor_date_broad <= last_encounter_date
      THEN incident_meningeal_tumor_date_broad
      ELSE last_encounter_date
    END AS followup_end_broad,

    DATE_DIFF(
      'day',
      index_date + INTERVAL %d DAY,
      CASE
        WHEN incident_meningioma_date_strict IS NOT NULL
         AND incident_meningioma_date_strict <= last_encounter_date
        THEN incident_meningioma_date_strict
        ELSE last_encounter_date
      END
    ) AS followup_days_strict,

    DATE_DIFF(
      'day',
      index_date + INTERVAL %d DAY,
      CASE
        WHEN incident_meningeal_tumor_date_broad IS NOT NULL
         AND incident_meningeal_tumor_date_broad <= last_encounter_date
        THEN incident_meningeal_tumor_date_broad
        ELSE last_encounter_date
      END
    ) AS followup_days_broad

FROM cohort_eligibility
WHERE %s;
",
OUTCOME_LAG_DAYS,
OUTCOME_LAG_DAYS,
OUTCOME_LAG_DAYS,
primary_where
))

try(dbExecute(con,
  "CREATE INDEX IF NOT EXISTS idx_analysis_primary_patient
   ON analysis_cohort_primary(patient_id);"
), silent = TRUE)

# ------------------------------------------------------------
# 11. QC: PRIMARY COHORT COUNTS + EVENTS
# ------------------------------------------------------------
primary_counts <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    COUNT(*) AS patients,

    SUM(meningioma_event_strict) AS strict_events,
    ROUND(100.0 * SUM(meningioma_event_strict) / COUNT(*), 4) AS strict_event_pct,

    SUM(meningeal_event_broad) AS broad_events,
    ROUND(100.0 * SUM(meningeal_event_broad) / COUNT(*), 4) AS broad_event_pct,

    MEDIAN(age_at_index) AS median_age,
    MEDIAN(followup_days_strict) / 365.25 AS median_followup_years,

    MIN(index_date) AS earliest_index,
    MAX(index_date) AS latest_index

FROM analysis_cohort_primary
GROUP BY exposure
ORDER BY exposure;
"))

fwrite(primary_counts, file.path(QC_DIR, "02_primary_counts_and_events.csv"))

# Index drug composition
drug_counts <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    index_drug,
    COUNT(*) AS patients,
    SUM(meningioma_event_strict) AS strict_events
FROM analysis_cohort_primary
GROUP BY exposure, index_drug
ORDER BY exposure, patients DESC;
"))

fwrite(drug_counts, file.path(QC_DIR, "02_primary_index_drugs.csv"))

# Calendar distribution
year_counts <- as.data.table(dbGetQuery(con, "
SELECT
    YEAR(index_date) AS index_year,
    exposure,
    COUNT(*) AS patients,
    SUM(meningioma_event_strict) AS strict_events
FROM analysis_cohort_primary
GROUP BY YEAR(index_date), exposure
ORDER BY index_year, exposure;
"))

fwrite(year_counts, file.path(QC_DIR, "02_primary_index_year.csv"))

# 5ARI strength distribution among primary cohort
strength_counts <- as.data.table(dbGetQuery(con, "
SELECT
    index_drug,
    COALESCE(index_finasteride_strength,
             index_dutasteride_strength,
             '[missing]') AS index_strength,
    COUNT(*) AS patients,
    SUM(meningioma_event_strict) AS strict_events
FROM analysis_cohort_primary
WHERE exposure='5ARI'
GROUP BY index_drug,
         COALESCE(index_finasteride_strength,
                  index_dutasteride_strength,
                  '[missing]')
ORDER BY index_drug, patients DESC;
"))

fwrite(strength_counts, file.path(QC_DIR, "02_primary_5ari_strengths.csv"))

# Baseline crude descriptive summary for early sanity checking.
baseline_crude <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    COUNT(*) AS n,
    AVG(age_at_index) AS mean_age,

    AVG(hypertension) AS hypertension_prop,
    AVG(diabetes) AS diabetes_prop,
    AVG(obesity) AS obesity_prop,
    AVG(ckd) AS ckd_prop,
    AVG(smoking) AS smoking_prop,
    AVG(cad) AS cad_prop,
    AVG(heart_failure) AS heart_failure_prop,
    AVG(cerebrovascular_disease) AS cerebrovascular_prop,
    AVG(hyperlipidemia) AS hyperlipidemia_prop,
    AVG(osteoporosis) AS osteoporosis_prop,
    AVG(alopecia) AS alopecia_prop,

    AVG(urinary_retention) AS urinary_retention_prop,
    AVG(uti) AS uti_prop,
    AVG(hematuria) AS hematuria_prop,
    AVG(bladder_stone) AS bladder_stone_prop,

    AVG(n_encounters_prior_365) AS mean_encounters_prior_365,
    AVG(n_ambulatory_prior_365) AS mean_ambulatory_prior_365,
    AVG(n_ed_prior_365) AS mean_ed_prior_365,
    AVG(n_inpatient_prior_365) AS mean_inpatient_prior_365

FROM analysis_cohort_primary
GROUP BY exposure
ORDER BY exposure;
"))

fwrite(baseline_crude, file.path(QC_DIR, "02_primary_baseline_crude.csv"))

# Follow-up quantiles
followup_qc <- as.data.table(dbGetQuery(con, "
SELECT
    exposure,
    QUANTILE_CONT(followup_days_strict / 365.25, 0.25) AS followup_q1_y,
    QUANTILE_CONT(followup_days_strict / 365.25, 0.50) AS followup_median_y,
    QUANTILE_CONT(followup_days_strict / 365.25, 0.75) AS followup_q3_y,
    QUANTILE_CONT(followup_days_strict / 365.25, 0.90) AS followup_p90_y,
    MAX(followup_days_strict / 365.25) AS followup_max_y
FROM analysis_cohort_primary
GROUP BY exposure
ORDER BY exposure;
"))

fwrite(followup_qc, file.path(QC_DIR, "02_primary_followup_qc.csv"))

# ------------------------------------------------------------
# 12. SENSITIVITY-READY FLAGS
# ------------------------------------------------------------
# Keep these in DuckDB now so later models do not require rescanning raw files.
dbExecute(con, "
CREATE OR REPLACE TABLE analysis_cohort_sensitivity AS
SELECT
    a.*,

    CASE
      WHEN exposure='5ARI'
       AND index_drug='finasteride'
       AND UPPER(TRIM(COALESCE(index_finasteride_strength,'')))='5 MG'
      THEN 1 ELSE 0
    END AS finasteride_5mg_only,

    CASE
      WHEN exposure='5ARI'
       AND index_drug='dutasteride'
       AND UPPER(TRIM(COALESCE(index_dutasteride_strength,'')))='0.5 MG'
      THEN 1 ELSE 0
    END AS dutasteride_05mg_only,

    CASE
      WHEN exposure='5ARI'
       AND index_drug='finasteride'
       AND UPPER(TRIM(COALESCE(index_finasteride_strength,'')))='1 MG'
      THEN 1 ELSE 0
    END AS finasteride_1mg,

    CASE
      WHEN followup_days_strict >= 365 THEN 1 ELSE 0
    END AS has_1y_post_landmark,

    CASE
      WHEN followup_days_strict >= 1095 THEN 1 ELSE 0
    END AS has_3y_post_landmark,

    CASE
      WHEN followup_days_strict >= 1826 THEN 1 ELSE 0
    END AS has_5y_post_landmark

FROM analysis_cohort_primary a;
")

# ------------------------------------------------------------
# 13. EXPORT COMPACT ANALYSIS TABLES
# ------------------------------------------------------------
primary_parquet <- sql_path(
  file.path(PARQUET_DIR, "analysis_cohort_primary.parquet")
)
sens_parquet <- sql_path(
  file.path(PARQUET_DIR, "analysis_cohort_sensitivity.parquet")
)

dbExecute(con, sprintf(
  "COPY analysis_cohort_primary TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000);",
  primary_parquet
))

dbExecute(con, sprintf(
  "COPY analysis_cohort_sensitivity TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000);",
  sens_parquet
))

# ------------------------------------------------------------
# 14. PRINT RESULTS
# ------------------------------------------------------------
cat("\n================ COHORT FLOW ======================\n")
print(flow)

cat("\n================ PRIMARY COUNTS / EVENTS ==========\n")
print(primary_counts)

cat("\n================ PRIMARY INDEX DRUGS ==============\n")
print(drug_counts)

cat("\n================ PRIMARY 5ARI STRENGTHS ===========\n")
print(strength_counts[1:min(.N, 30)])

cat("\n================ FOLLOW-UP ========================\n")
print(followup_qc)

cat("\n================ CRUDE BASELINE ===================\n")
print(baseline_crude)

cat("\nStage 2 complete.\n")
cat("QC directory:", QC_DIR, "\n")
cat("Primary cohort:", primary_parquet, "\n")
cat("Sensitivity cohort:", sens_parquet, "\n")
cat("DuckDB:", DB_PATH, "\n")

writeLines(
  c(
    paste0("INDEX_START=", INDEX_START),
    paste0("DATA_CUTOFF=", data_cutoff),
    paste0("INDEX_END=", INDEX_END),
    paste0("MIN_AGE=", MIN_AGE),
    paste0("MAX_AGE=", MAX_AGE),
    paste0("LOOKBACK_DAYS=", LOOKBACK_DAYS),
    paste0("BASELINE_COVARIATE_DAYS=", BASELINE_COVARIATE_DAYS),
    paste0("BPH_LOOKBACK_DAYS=", BPH_LOOKBACK_DAYS),
    paste0("BPH_GRACE_DAYS=", BPH_GRACE_DAYS),
    paste0("TREATMENT_GRACE_DAYS=", TREATMENT_GRACE_DAYS),
    paste0("OUTCOME_LAG_DAYS=", OUTCOME_LAG_DAYS),
    paste0("MIN_POST_LAG_DAYS=", MIN_POST_LAG_DAYS),
    paste0("EXCLUDE_KNOWN_FINASTERIDE_1MG=", EXCLUDE_KNOWN_FINASTERIDE_1MG)
  ),
  file.path(QC_DIR, "02_design_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(QC_DIR, "02_sessionInfo.txt")
)

dbDisconnect(con, shutdown = TRUE)
