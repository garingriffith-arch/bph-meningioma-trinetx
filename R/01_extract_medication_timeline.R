# ============================================================
# 01_extract_medication_timeline.R
# Finasteride / Dutasteride and Incident Meningioma Study
#
# Purpose:
#   1) Scan medication_ingredient.csv ONCE.
#   2) Keep only the seven medications from the source query.
#   3) Deduplicate study-medication events.
#   4) Build patient-level first exposure dates.
#   5) Define preliminary clean index candidates.
#   6) Produce QC counts, calendar trends, strengths, and crossover data.
#
# Run AFTER 00_validate_export.R.
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

required_pkgs <- c("DBI", "duckdb", "data.table")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

library(DBI)
library(duckdb)
library(data.table)

# -----------------------------
# 1. Paths / performance
# -----------------------------
RAW_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/Data/BPH"
R_DIR   <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"

WORK_DIR <- "C:/TriNetX/bph_m_work"
DB_DIR   <- file.path(WORK_DIR, "duckdb")
TEMP_DIR <- file.path(WORK_DIR, "duckdb_temp")
QC_DIR   <- file.path(WORK_DIR, "qc")
PARQUET_DIR <- file.path(WORK_DIR, "parquet")

dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PARQUET_DIR, recursive = TRUE, showWarnings = FALSE)

DB_PATH <- file.path(DB_DIR, "bph_m.duckdb")
MED_INGREDIENT <- file.path(RAW_DIR, "medication_ingredient.csv")

if (!file.exists(MED_INGREDIENT)) {
  stop("medication_ingredient.csv not found at: ", MED_INGREDIENT)
}

N_THREADS <- min(24L, max(4L, parallel::detectCores(logical = TRUE) - 6L))
MEMORY_LIMIT <- "48GB"

sql_path <- function(x) {
  x <- normalizePath(x, winslash = "/", mustWork = FALSE)
  gsub("'", "''", x, fixed = TRUE)
}

# -----------------------------
# 2. Study medication dictionary
# -----------------------------
# Exact RxNorm ingredient codes from your TriNetX query.
drug_dictionary <- data.table(
  code = c(
    "25025",   # finasteride
    "228790",  # dutasteride
    "77492",   # tamsulosin
    "17300",   # alfuzosin
    "720825",  # silodosin
    "49276",   # doxazosin
    "37798"    # terazosin
  ),
  drug = c(
    "finasteride",
    "dutasteride",
    "tamsulosin",
    "alfuzosin",
    "silodosin",
    "doxazosin",
    "terazosin"
  ),
  class = c(
    "5ARI",
    "5ARI",
    "URO_ALPHA",
    "URO_ALPHA",
    "URO_ALPHA",
    "NONURO_ALPHA",
    "NONURO_ALPHA"
  )
)

fwrite(drug_dictionary, file.path(QC_DIR, "01_drug_dictionary.csv"))

code_sql <- paste(sprintf("'%s'", drug_dictionary$code), collapse = ", ")

# -----------------------------
# 3. Connect to DuckDB
# -----------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = DB_PATH, read_only = FALSE)

on.exit({
  try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

dbExecute(con, sprintf("PRAGMA threads=%d;", N_THREADS))
dbExecute(con, sprintf("PRAGMA memory_limit='%s';", MEMORY_LIMIT))
dbExecute(con, sprintf("PRAGMA temp_directory='%s';", sql_path(TEMP_DIR)))
dbExecute(con, "SET preserve_insertion_order=false;")

# Confirm script 00 ran.
needed_tables <- c("patient_core", "patient_cohort_core")
existing <- dbListTables(con)

if (!all(needed_tables %in% existing)) {
  stop(
    "Run 00_validate_export.R first. Missing DuckDB tables: ",
    paste(setdiff(needed_tables, existing), collapse = ", ")
  )
}

# -----------------------------
# 4. Scan the 229-GB medication ingredient file ONCE
# -----------------------------
# We retain route/brand/strength because:
#   - finasteride 1 mg is commonly used for alopecia;
#   - finasteride 5 mg is the usual BPH dose;
#   - this lets us later build dose/indication sensitivity analyses.
#
# start_date is converted from TriNetX YYYYMMDD text to a true DATE.
#
# source_id is retained so we can later inspect EHR vs Closed Claims A.
cat("\nStarting single-pass medication extraction.\n")
cat("This is the heavy step. Do not interrupt it unless necessary.\n")
cat("Source:", MED_INGREDIENT, "\n\n")

t0 <- Sys.time()

dbExecute(con, sprintf(
"
CREATE OR REPLACE TABLE study_med_events_raw AS
SELECT
    patient_id,
    encounter_id,
    unique_id,
    code AS drug_code,
    CASE code
        WHEN '25025'  THEN 'finasteride'
        WHEN '228790' THEN 'dutasteride'
        WHEN '77492'  THEN 'tamsulosin'
        WHEN '17300'  THEN 'alfuzosin'
        WHEN '720825' THEN 'silodosin'
        WHEN '49276'  THEN 'doxazosin'
        WHEN '37798'  THEN 'terazosin'
    END AS drug,
    CASE
        WHEN code IN ('25025','228790') THEN '5ARI'
        WHEN code IN ('77492','17300','720825') THEN 'URO_ALPHA'
        WHEN code IN ('49276','37798') THEN 'NONURO_ALPHA'
    END AS drug_class,
    TRY_CAST(STRPTIME(start_date, '%%Y%%m%%d') AS DATE) AS start_date,
    route,
    brand,
    strength,
    source_id
FROM read_csv_auto(
    '%s',
    header=true,
    all_varchar=true,
    ignore_errors=false
)
WHERE code_system = 'RxNorm'
  AND code IN (%s);
",
sql_path(MED_INGREDIENT),
code_sql
))

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat("\nMedication raw extraction completed in",
    round(as.numeric(elapsed), 1), "minutes.\n")

# -----------------------------
# 5. QC invalid/missing dates
# -----------------------------
date_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    COUNT(*) AS rows_total,
    SUM(CASE WHEN start_date IS NULL THEN 1 ELSE 0 END) AS missing_or_invalid_start_date
FROM study_med_events_raw;
"))

fwrite(date_qc, file.path(QC_DIR, "01_medication_date_qc.csv"))
print(date_qc)

# -----------------------------
# 6. Deduplicate medication events
# -----------------------------
# Multiple source records can represent the same observed medication event.
# For exposure timing, one patient/drug/date is sufficient.
#
# We preserve a representative route/brand/strength/source value using MIN().
dbExecute(con,
"
CREATE OR REPLACE TABLE study_med_events AS
SELECT
    patient_id,
    drug_code,
    drug,
    drug_class,
    start_date,
    MIN(route) AS route,
    MIN(brand) AS brand,
    MIN(strength) AS strength,
    MIN(source_id) AS source_id,
    COUNT(*) AS raw_rows_collapsed
FROM study_med_events_raw
WHERE start_date IS NOT NULL
GROUP BY
    patient_id,
    drug_code,
    drug,
    drug_class,
    start_date;
")

# Indexes improve later patient/date joins inside DuckDB.
try(dbExecute(con,
  "CREATE INDEX IF NOT EXISTS idx_study_med_patient ON study_med_events(patient_id);"
), silent = TRUE)

try(dbExecute(con,
  "CREATE INDEX IF NOT EXISTS idx_study_med_patient_date ON study_med_events(patient_id, start_date);"
), silent = TRUE)

# -----------------------------
# 7. Patient-level first exposure dates
# -----------------------------
dbExecute(con,
"
CREATE OR REPLACE TABLE study_med_firsts AS
SELECT
    patient_id,

    MIN(start_date) AS first_any_study_med,

    MIN(CASE WHEN drug = 'finasteride' THEN start_date END) AS first_finasteride,
    MIN(CASE WHEN drug = 'dutasteride' THEN start_date END) AS first_dutasteride,

    MIN(CASE WHEN drug = 'tamsulosin' THEN start_date END) AS first_tamsulosin,
    MIN(CASE WHEN drug = 'alfuzosin' THEN start_date END) AS first_alfuzosin,
    MIN(CASE WHEN drug = 'silodosin' THEN start_date END) AS first_silodosin,
    MIN(CASE WHEN drug = 'doxazosin' THEN start_date END) AS first_doxazosin,
    MIN(CASE WHEN drug = 'terazosin' THEN start_date END) AS first_terazosin,

    MIN(CASE WHEN drug_class = '5ARI' THEN start_date END) AS first_5ari,
    MIN(CASE WHEN drug_class = 'URO_ALPHA' THEN start_date END) AS first_uro_alpha,
    MIN(CASE WHEN drug_class = 'NONURO_ALPHA' THEN start_date END) AS first_nonuro_alpha,

    COUNT(*) AS n_unique_med_dates,
    COUNT(DISTINCT drug_code) AS n_distinct_study_drugs

FROM study_med_events
GROUP BY patient_id;
")

# -----------------------------
# 8. Identify what was started on the FIRST observed study-med date
# -----------------------------
# This intentionally treats same-day multi-class starts as combination/ambiguous.
dbExecute(con,
"
CREATE OR REPLACE TABLE first_day_classes AS
SELECT
    f.patient_id,
    f.first_any_study_med AS index_date_prelim,

    COUNT(DISTINCT e.drug_class) AS n_classes_on_first_day,
    COUNT(DISTINCT e.drug) AS n_drugs_on_first_day,

    MAX(CASE WHEN e.drug_class = '5ARI' THEN 1 ELSE 0 END) AS first_day_5ari,
    MAX(CASE WHEN e.drug_class = 'URO_ALPHA' THEN 1 ELSE 0 END) AS first_day_uro_alpha,
    MAX(CASE WHEN e.drug_class = 'NONURO_ALPHA' THEN 1 ELSE 0 END) AS first_day_nonuro_alpha,

    STRING_AGG(DISTINCT e.drug, ', ' ORDER BY e.drug) AS first_day_drugs

FROM study_med_firsts f
JOIN study_med_events e
  ON e.patient_id = f.patient_id
 AND e.start_date = f.first_any_study_med

GROUP BY
    f.patient_id,
    f.first_any_study_med;
")

# -----------------------------
# 9. Preliminary exposure candidates
# -----------------------------
# PRIMARY candidate logic:
#   5ARI:
#      first observed study-medication day has ONLY 5ARI medication(s)
#   URO_ALPHA:
#      first observed study-medication day has ONLY uroselective alpha blocker(s)
#
# NONURO alpha starts and combination starts are retained but not eligible for
# the eventual primary active-comparator cohort.
#
# IMPORTANT:
# This is still PRELIMINARY. We have NOT yet imposed:
#   - >=365 days prior observable history
#   - BPH diagnosis timing relative to index
#   - prior meningioma exclusion
#   - prostate cancer / androgen deprivation exclusion
#   - age >=40 at index
# Those come in the next cohort-construction script.
dbExecute(con,
"
CREATE OR REPLACE TABLE index_candidates_prelim AS
SELECT
    f.*,
    d.index_date_prelim,
    d.n_classes_on_first_day,
    d.n_drugs_on_first_day,
    d.first_day_drugs,

    CASE
        WHEN d.n_classes_on_first_day = 1 AND d.first_day_5ari = 1
            THEN '5ARI'
        WHEN d.n_classes_on_first_day = 1 AND d.first_day_uro_alpha = 1
            THEN 'URO_ALPHA'
        WHEN d.n_classes_on_first_day = 1 AND d.first_day_nonuro_alpha = 1
            THEN 'NONURO_ALPHA'
        ELSE 'COMBINATION_OR_AMBIGUOUS'
    END AS exposure_prelim

FROM study_med_firsts f
JOIN first_day_classes d USING (patient_id);
")

# -----------------------------
# 10. Crossover / treatment sequencing
# -----------------------------
# Days between class starts. These will later support:
#   - treatment switching description
#   - excluding near-simultaneous starts
#   - per-protocol / as-treated sensitivities
dbExecute(con,
"
CREATE OR REPLACE TABLE class_sequence AS
SELECT
    patient_id,
    first_5ari,
    first_uro_alpha,
    first_nonuro_alpha,

    CASE
        WHEN first_5ari IS NOT NULL AND first_uro_alpha IS NOT NULL
        THEN DATE_DIFF('day', first_5ari, first_uro_alpha)
    END AS days_5ari_to_uro_alpha,

    CASE
        WHEN first_uro_alpha IS NOT NULL AND first_5ari IS NOT NULL
        THEN DATE_DIFF('day', first_uro_alpha, first_5ari)
    END AS days_uro_alpha_to_5ari

FROM study_med_firsts;
")

# -----------------------------
# 11. QC: counts by medication
# -----------------------------
med_counts <- as.data.table(dbGetQuery(con,
"
SELECT
    drug,
    drug_class,
    COUNT(*) AS unique_patient_drug_dates,
    COUNT(DISTINCT patient_id) AS patients,
    MIN(start_date) AS earliest_date,
    MAX(start_date) AS latest_date
FROM study_med_events
GROUP BY drug, drug_class
ORDER BY patients DESC;
"))

fwrite(med_counts, file.path(QC_DIR, "01_medication_counts.csv"))

# -----------------------------
# 12. QC: preliminary exposure counts
# -----------------------------
exposure_counts <- as.data.table(dbGetQuery(con,
"
SELECT
    exposure_prelim,
    COUNT(*) AS patients,
    MIN(index_date_prelim) AS earliest_index,
    MAX(index_date_prelim) AS latest_index
FROM index_candidates_prelim
GROUP BY exposure_prelim
ORDER BY patients DESC;
"))

fwrite(exposure_counts, file.path(QC_DIR, "01_preliminary_exposure_counts.csv"))

# -----------------------------
# 13. QC: calendar year
# -----------------------------
calendar_counts <- as.data.table(dbGetQuery(con,
"
SELECT
    YEAR(index_date_prelim) AS index_year,
    exposure_prelim,
    COUNT(*) AS patients
FROM index_candidates_prelim
WHERE index_date_prelim IS NOT NULL
GROUP BY YEAR(index_date_prelim), exposure_prelim
ORDER BY index_year, exposure_prelim;
"))

fwrite(calendar_counts, file.path(QC_DIR, "01_index_year_counts.csv"))

# -----------------------------
# 14. QC: first-day medication combinations
# -----------------------------
first_day_combo_counts <- as.data.table(dbGetQuery(con,
"
SELECT
    exposure_prelim,
    first_day_drugs,
    COUNT(*) AS patients
FROM index_candidates_prelim
GROUP BY exposure_prelim, first_day_drugs
ORDER BY patients DESC;
"))

fwrite(first_day_combo_counts, file.path(QC_DIR, "01_first_day_combinations.csv"))

# -----------------------------
# 15. QC: finasteride / dutasteride strength
# -----------------------------
# This is crucial for assessing whether we can distinguish BPH-dose finasteride
# from alopecia-dose finasteride.
strength_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    drug,
    COALESCE(NULLIF(TRIM(strength), ''), '[missing]') AS strength,
    COUNT(*) AS unique_patient_drug_dates,
    COUNT(DISTINCT patient_id) AS patients
FROM study_med_events
WHERE drug IN ('finasteride', 'dutasteride')
GROUP BY drug, COALESCE(NULLIF(TRIM(strength), ''), '[missing]')
ORDER BY drug, patients DESC;
"))

fwrite(strength_qc, file.path(QC_DIR, "01_5ari_strengths.csv"))

# -----------------------------
# 16. QC: source type
# -----------------------------
source_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    drug,
    source_id,
    COUNT(*) AS unique_patient_drug_dates,
    COUNT(DISTINCT patient_id) AS patients
FROM study_med_events
GROUP BY drug, source_id
ORDER BY drug, patients DESC;
"))

fwrite(source_qc, file.path(QC_DIR, "01_medication_source_qc.csv"))

# -----------------------------
# 17. QC: overlap/crossover
# -----------------------------
sequence_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    COUNT(*) AS patients_with_any_study_med,

    SUM(CASE WHEN first_5ari IS NOT NULL THEN 1 ELSE 0 END) AS ever_5ari,
    SUM(CASE WHEN first_uro_alpha IS NOT NULL THEN 1 ELSE 0 END) AS ever_uro_alpha,
    SUM(CASE WHEN first_nonuro_alpha IS NOT NULL THEN 1 ELSE 0 END) AS ever_nonuro_alpha,

    SUM(CASE WHEN first_5ari IS NOT NULL AND first_uro_alpha IS NOT NULL THEN 1 ELSE 0 END)
        AS ever_both_5ari_and_uro_alpha,

    SUM(CASE
        WHEN first_5ari IS NOT NULL
         AND first_uro_alpha IS NOT NULL
         AND ABS(DATE_DIFF('day', first_5ari, first_uro_alpha)) <= 30
        THEN 1 ELSE 0 END) AS class_starts_within_30d,

    SUM(CASE
        WHEN first_5ari IS NOT NULL
         AND first_uro_alpha IS NOT NULL
         AND ABS(DATE_DIFF('day', first_5ari, first_uro_alpha)) <= 90
        THEN 1 ELSE 0 END) AS class_starts_within_90d

FROM study_med_firsts;
"))

fwrite(sequence_qc, file.path(QC_DIR, "01_class_sequence_qc.csv"))

# -----------------------------
# 18. Export compact derived tables to Parquet
# -----------------------------
# These are dramatically smaller/faster than rescanning 229 GB of CSV.
parquet_events <- sql_path(file.path(PARQUET_DIR, "study_med_events.parquet"))
parquet_firsts <- sql_path(file.path(PARQUET_DIR, "study_med_firsts.parquet"))
parquet_index  <- sql_path(file.path(PARQUET_DIR, "index_candidates_prelim.parquet"))

dbExecute(con, sprintf(
  "COPY study_med_events TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000);",
  parquet_events
))

dbExecute(con, sprintf(
  "COPY study_med_firsts TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000);",
  parquet_firsts
))

dbExecute(con, sprintf(
  "COPY index_candidates_prelim TO '%s'
   (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000);",
  parquet_index
))

# -----------------------------
# 19. Print the high-value QC outputs
# -----------------------------
cat("\n================ MEDICATION COUNTS ================\n")
print(med_counts)

cat("\n================ PRELIM EXPOSURE COUNTS ============\n")
print(exposure_counts)

cat("\n================ CLASS SEQUENCE QC =================\n")
print(sequence_qc)

cat("\n================ 5ARI STRENGTHS (TOP 30) ===========\n")
print(strength_qc[1:min(.N, 30)])

cat("\n================ FIRST-DAY COMBINATIONS (TOP 30) ===\n")
print(first_day_combo_counts[1:min(.N, 30)])

cat("\nMedication timeline extraction complete.\n")
cat("QC output:", QC_DIR, "\n")
cat("Parquet output:", PARQUET_DIR, "\n")
cat("DuckDB:", DB_PATH, "\n")

# Save script-level session info.
writeLines(
  capture.output(sessionInfo()),
  file.path(QC_DIR, "01_sessionInfo.txt")
)

dbDisconnect(con, shutdown = TRUE)
