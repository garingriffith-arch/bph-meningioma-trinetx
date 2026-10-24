# ============================================================
# 00_validate_export.R
# Finasteride / Dutasteride and Incident Meningioma Study
# TriNetX Dataworks export validation + local DuckDB setup
#
# Run this FIRST.
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

# -----------------------------
# 0. Packages
# -----------------------------
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
# 1. Paths
# -----------------------------
RAW_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/Data/BPH"
R_DIR   <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"

# Keep the database / spill files OUTSIDE OneDrive.
# C: is your NVMe drive from the machine screenshot.
WORK_DIR <- "C:/TriNetX/bph_m_work"
DB_DIR   <- file.path(WORK_DIR, "duckdb")
TEMP_DIR <- file.path(WORK_DIR, "duckdb_temp")
QC_DIR   <- file.path(WORK_DIR, "qc")
PARQUET_DIR <- file.path(WORK_DIR, "parquet")

dir.create(R_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PARQUET_DIR, recursive = TRUE, showWarnings = FALSE)

DB_PATH <- file.path(DB_DIR, "bph_m.duckdb")

# -----------------------------
# 2. Performance settings
# -----------------------------
# i9-14900KF: 32 logical threads. Leave several for Windows/OneDrive/RStudio.
N_THREADS <- min(24L, max(4L, parallel::detectCores(logical = TRUE) - 6L))

# 64 GB RAM system: leave ~16 GB for Windows + RStudio + filesystem cache.
MEMORY_LIMIT <- "48GB"

sql_path <- function(x) {
  x <- normalizePath(x, winslash = "/", mustWork = FALSE)
  gsub("'", "''", x, fixed = TRUE)
}

qident <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))

# -----------------------------
# 3. Check expected files
# -----------------------------
expected_files <- c(
  "chemo_lines.csv",
  "cohort_details.csv",
  "datadictionary.pdf",
  "datadictionary.xlsx",
  "dataset_details.csv",
  "diagnosis.csv",
  "encounter.csv",
  "genomic.csv",
  "lab_result.csv",
  "manifest.csv",
  "medication_drug.csv",
  "medication_ingredient.csv",
  "oncology_treatment.csv",
  "patient.csv",
  "patient_cohort.csv",
  "procedure.csv",
  "standardized_terminology.csv",
  "tumor.csv",
  "tumor_properties.csv",
  "vitals_signs.csv"
)

present_names <- list.files(RAW_DIR, full.names = FALSE)
file_check <- data.table(
  file = expected_files,
  exists = expected_files %in% present_names
)

# Accept datadictionary files with suffixes/renamed copies; they are not needed by code.
file_check[file %in% c("datadictionary.pdf", "datadictionary.xlsx"), exists := TRUE]

if (any(!file_check$exists)) {
  print(file_check[exists == FALSE])
  stop("One or more required TriNetX CSV files are missing from RAW_DIR.")
}

# Logical file sizes. OneDrive placeholders can report size even if not hydrated,
# so this is a QC summary, not proof that every byte is resident locally.
csv_paths <- list.files(RAW_DIR, pattern = "\\.csv$", full.names = TRUE)
sizes <- file.info(csv_paths)$size

file_size_qc <- data.table(
  file = basename(csv_paths),
  bytes = as.numeric(sizes),
  GiB = round(as.numeric(sizes) / 1024^3, 3)
)[order(-GiB)]

fwrite(file_size_qc, file.path(QC_DIR, "00_file_sizes.csv"))
cat("\nLargest files (GiB):\n")
print(file_size_qc[1:min(.N, 20)])

cat("\nTotal CSV logical size:",
    round(sum(file_size_qc$bytes, na.rm = TRUE) / 1024^3, 1),
    "GiB\n")

# -----------------------------
# 4. Open DuckDB
# -----------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = DB_PATH, read_only = FALSE)

on.exit({
  try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

dbExecute(con, sprintf("PRAGMA threads=%d;", N_THREADS))
dbExecute(con, sprintf("PRAGMA memory_limit='%s';", MEMORY_LIMIT))
dbExecute(con, sprintf("PRAGMA temp_directory='%s';", sql_path(TEMP_DIR)))
dbExecute(con, "SET preserve_insertion_order=false;")

cat("\nDuckDB database:", DB_PATH, "\n")
cat("DuckDB threads:", N_THREADS, "\n")
cat("DuckDB memory limit:", MEMORY_LIMIT, "\n")
cat("DuckDB temp:", TEMP_DIR, "\n")

# -----------------------------
# 5. Helper: read TriNetX CSV as all-character
# -----------------------------
csv_scan <- function(filename) {
  p <- sql_path(file.path(RAW_DIR, filename))
  sprintf(
    "read_csv_auto('%s', header=true, all_varchar=true, ignore_errors=false)",
    p
  )
}

# -----------------------------
# 6. Load the SMALL metadata tables
# -----------------------------
small_tables <- c(
  cohort_details = "cohort_details.csv",
  dataset_details = "dataset_details.csv",
  manifest = "manifest.csv"
)

for (tbl in names(small_tables)) {
  sql <- sprintf(
    "CREATE OR REPLACE TABLE %s AS SELECT * FROM %s;",
    tbl, csv_scan(small_tables[[tbl]])
  )
  dbExecute(con, sql)
}

# Standardized terminology is relatively small (~hundreds of MB) and useful
# throughout the project. Materialize it once locally.
dbExecute(con, sprintf(
  "CREATE OR REPLACE TABLE standardized_terminology AS
   SELECT * FROM %s;",
  csv_scan("standardized_terminology.csv")
))

# -----------------------------
# 7. Materialize patient table
# -----------------------------
# The Dataworks dictionary defines patient_id as the patient-table primary key.
# We explicitly cast year_of_birth and retain source/geography fields.
dbExecute(con, sprintf(
"
CREATE OR REPLACE TABLE patient_core AS
SELECT
    patient_id,
    sex,
    race,
    ethnicity,
    marital_status,
    TRY_CAST(year_of_birth AS INTEGER) AS year_of_birth,
    reason_yob_missing,
    TRY_CAST(month_year_death AS BIGINT) AS month_year_death,
    death_date_source_id,
    patient_regional_location,
    source_id
FROM %s;
",
csv_scan("patient.csv")
))

# Patient-cohort membership can be useful for reconciling query cohort counts.
dbExecute(con, sprintf(
"
CREATE OR REPLACE TABLE patient_cohort_core AS
SELECT
    patient_id,
    cohort_name,
    TRY_CAST(cohort_number AS BIGINT) AS cohort_number,
    source_id
FROM %s;
",
csv_scan("patient_cohort.csv")
))

# -----------------------------
# 8. QC summaries
# -----------------------------
dataset_details <- as.data.table(dbGetQuery(con, "SELECT * FROM dataset_details;"))
cohort_details  <- as.data.table(dbGetQuery(con, "SELECT * FROM cohort_details;"))
manifest        <- as.data.table(dbGetQuery(con, "SELECT * FROM manifest;"))

fwrite(dataset_details, file.path(QC_DIR, "00_dataset_details.csv"))
fwrite(cohort_details,  file.path(QC_DIR, "00_cohort_details.csv"))
fwrite(manifest,        file.path(QC_DIR, "00_manifest.csv"))

patient_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    COUNT(*) AS patient_rows,
    COUNT(DISTINCT patient_id) AS distinct_patient_ids,
    COUNT(*) - COUNT(DISTINCT patient_id) AS duplicate_patient_rows,
    SUM(CASE WHEN sex = 'M' THEN 1 ELSE 0 END) AS male_rows,
    SUM(CASE WHEN sex = 'F' THEN 1 ELSE 0 END) AS female_rows,
    SUM(CASE WHEN sex IS NULL OR sex = '' OR sex = 'Unknown' THEN 1 ELSE 0 END) AS unknown_sex_rows,
    SUM(CASE WHEN year_of_birth IS NULL THEN 1 ELSE 0 END) AS missing_yob_rows,
    MIN(year_of_birth) AS min_yob,
    MAX(year_of_birth) AS max_yob
FROM patient_core;
"))

fwrite(patient_qc, file.path(QC_DIR, "00_patient_qc.csv"))

patient_cohort_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    cohort_number,
    cohort_name,
    COUNT(*) AS rows,
    COUNT(DISTINCT patient_id) AS distinct_patient_ids
FROM patient_cohort_core
GROUP BY cohort_number, cohort_name
ORDER BY cohort_number, cohort_name;
"))

fwrite(patient_cohort_qc, file.path(QC_DIR, "00_patient_cohort_qc.csv"))

source_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    source_id,
    COUNT(*) AS n
FROM patient_core
GROUP BY source_id
ORDER BY n DESC;
"))

fwrite(source_qc, file.path(QC_DIR, "00_patient_source_qc.csv"))

race_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    race,
    COUNT(*) AS n
FROM patient_core
GROUP BY race
ORDER BY n DESC;
"))

fwrite(race_qc, file.path(QC_DIR, "00_race_qc.csv"))

ethnicity_qc <- as.data.table(dbGetQuery(con,
"
SELECT
    ethnicity,
    COUNT(*) AS n
FROM patient_core
GROUP BY ethnicity
ORDER BY n DESC;
"))

fwrite(ethnicity_qc, file.path(QC_DIR, "00_ethnicity_qc.csv"))

# -----------------------------
# 9. Basic integrity checks
# -----------------------------
cat("\n================ DATASET DETAILS ================\n")
print(dataset_details)

cat("\n================ COHORT DETAILS =================\n")
print(cohort_details)

cat("\n================ PATIENT QC =====================\n")
print(patient_qc)

cat("\n================ PATIENT COHORT QC ==============\n")
print(patient_cohort_qc)

cat("\n================ LARGEST RAW FILES ==============\n")
print(file_size_qc[1:min(.N, 12)])

if (patient_qc$duplicate_patient_rows[1] != 0) {
  warning("patient.csv contains duplicate patient_id rows; inspect before proceeding.")
}

if (patient_qc$female_rows[1] > 0) {
  warning("Female rows are present despite the male Query Builder filter. Inspect query/export.")
}

# -----------------------------
# 10. Save reproducibility info
# -----------------------------
session_txt <- capture.output(sessionInfo())
writeLines(session_txt, file.path(QC_DIR, "00_sessionInfo.txt"))

config_txt <- c(
  paste0("RAW_DIR=", RAW_DIR),
  paste0("R_DIR=", R_DIR),
  paste0("WORK_DIR=", WORK_DIR),
  paste0("DB_PATH=", DB_PATH),
  paste0("TEMP_DIR=", TEMP_DIR),
  paste0("N_THREADS=", N_THREADS),
  paste0("MEMORY_LIMIT=", MEMORY_LIMIT)
)
writeLines(config_txt, file.path(QC_DIR, "00_runtime_config.txt"))

cat("\nValidation complete.\n")
cat("QC output:", QC_DIR, "\n")
cat("DuckDB:", DB_PATH, "\n")
cat("\nNEXT: run 01_extract_medication_timeline.R\n")

dbDisconnect(con, shutdown = TRUE)
