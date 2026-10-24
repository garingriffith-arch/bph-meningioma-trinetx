# NOTE: This stage is wrapped in a function so on.exit() is tied to
# the stage function rather than source()/eval(). This prevents the
# DuckDB connection from being closed immediately when sourced.

run_stage_08 <- function() {
  # ============================================================
  # 08_build_permissive_and_latency_cohorts.R
  #
  # Builds the deliberately permissive treatment-strategy sensitivity
  # requested for the finasteride / meningioma project.
  #
  # Exposure strategy:
  #   5ARI:
  #     ANY first-day regimen containing finasteride or dutasteride,
  #     including same-day alpha-blocker combination therapy.
  #
  #   URO_ALPHA:
  #     first-day regimen containing a uroselective alpha blocker
  #     and NO 5ARI. Same-day nonuroselective alpha use is allowed.
  #
  # Two permissive cohorts:
  #   A) MAXIMAL DEFENSIBLE:
  #      does NOT require BPH near index, does NOT exclude prostate cancer,
  #      allows 1-mg/unknown finasteride, combination therapy, later crossover.
  #
  #   B) PERMISSIVE BPH:
  #      same broad exposure strategy, but requires BPH near index and
  #      excludes prior prostate cancer.
  #
  # Also builds a CLEAN zero-lag latency cohort so that 0-1, 1-2, 2-5,
  # and >=5-year interval-specific estimates can be compared with Giraldi.
  #
  # Heavy scans:
  #   encounter.csv once for only patients not already summarized
  #   diagnosis.csv once for the whole permissive candidate set
  #
  # Both are checkpointed.
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

  ENCOUNTER_FILE <- file.path(RAW_DIR, "encounter.csv")
  DIAGNOSIS_FILE <- file.path(RAW_DIR, "diagnosis.csv")

  dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(PARQUET_DIR, recursive = TRUE, showWarnings = FALSE)

  N_THREADS <- min(24L, max(4L, parallel::detectCores(logical = TRUE) - 6L))
  MEMORY_LIMIT <- "48GB"
  REUSE_EXISTING_HEAVY_TABLES <- TRUE

  MIN_AGE <- 40L
  MAX_AGE <- 89L
  INDEX_START <- as.Date("2004-01-01")
  LOOKBACK_DAYS <- 365L
  BASELINE_DAYS <- 1095L
  BPH_LOOKBACK_DAYS <- 365L
  BPH_GRACE_DAYS <- 30L

  sql_path <- function(x) {
    x <- normalizePath(x, winslash = "/", mustWork = FALSE)
    gsub("'", "''", x, fixed = TRUE)
  }

  sql_fill <- function(template, values) {
    out <- template
    for (nm in names(values)) {
      out <- gsub(
        paste0("{{", nm, "}}"),
        as.character(values[[nm]]),
        out,
        fixed = TRUE
      )
    }
    out
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
    "patient_core",
    "study_med_firsts",
    "study_med_events",
    "first_day_classes",
    "candidate_observation"
  )
  missing <- setdiff(needed, dbListTables(con))
  if (length(missing) > 0L) {
    stop("Missing required prior-stage tables: ", paste(missing, collapse=", "))
  }

  # ------------------------------------------------------------
  # 1. DATA CUTOFF
  # ------------------------------------------------------------
  data_cutoff <- as.Date(
    dbGetQuery(
      con,
      "SELECT MAX(last_encounter_date) AS x FROM candidate_observation;"
    )$x[1]
  )

  if (is.na(data_cutoff)) stop("Could not derive data cutoff.")

  INDEX_END_LAG1 <- data_cutoff - 365L
  INDEX_END_ZERO <- data_cutoff - 30L

  cat("\nObserved data cutoff:", format(data_cutoff), "\n")
  cat("Permissive 1-y-lag index end:", format(INDEX_END_LAG1), "\n")
  cat("Zero-lag latency index end:", format(INDEX_END_ZERO), "\n")

  # ------------------------------------------------------------
  # 2. PERMISSIVE INDEX STRATEGIES
  # ------------------------------------------------------------
  dbExecute(con, "
  CREATE OR REPLACE TABLE permissive_index0 AS
  SELECT
    f.patient_id,
    f.first_any_study_med AS index_date,

    CASE
      WHEN d.first_day_5ari=1
        THEN '5ARI'
      WHEN d.first_day_5ari=0
       AND d.first_day_uro_alpha=1
        THEN 'URO_ALPHA'
      ELSE NULL
    END AS exposure,

    d.first_day_drugs AS index_drug,
    d.n_classes_on_first_day,
    d.n_drugs_on_first_day,
    d.first_day_5ari,
    d.first_day_uro_alpha,
    d.first_day_nonuro_alpha,

    f.first_finasteride,
    f.first_dutasteride,
    f.first_tamsulosin,
    f.first_alfuzosin,
    f.first_silodosin,
    f.first_doxazosin,
    f.first_terazosin,
    f.first_5ari,
    f.first_uro_alpha,
    f.first_nonuro_alpha,

    p.sex,
    p.race,
    p.ethnicity,
    p.marital_status,
    p.year_of_birth,
    p.patient_regional_location,
    p.source_id AS patient_source_id,

    YEAR(f.first_any_study_med) - p.year_of_birth AS age_at_index

  FROM study_med_firsts f
  JOIN first_day_classes d USING (patient_id)
  JOIN patient_core p USING (patient_id)

  WHERE
    (
      d.first_day_5ari=1
      OR (
        d.first_day_5ari=0
        AND d.first_day_uro_alpha=1
      )
    );
  ")

  # Index-day finasteride strength / 1-mg flag.
  dbExecute(con, "
  CREATE OR REPLACE TABLE permissive_index_strength AS
  SELECT
    p.patient_id,

    MIN(
      CASE
        WHEN e.drug='finasteride'
         AND e.start_date=p.index_date
        THEN e.strength
      END
    ) AS index_finasteride_strength,

    MIN(
      CASE
        WHEN e.drug='dutasteride'
         AND e.start_date=p.index_date
        THEN e.strength
      END
    ) AS index_dutasteride_strength

  FROM permissive_index0 p
  LEFT JOIN study_med_events e
    ON e.patient_id=p.patient_id
   AND e.start_date=p.index_date

  GROUP BY p.patient_id;
  ")

  # ------------------------------------------------------------
  # 3. OBSERVATION / UTILIZATION
  # ------------------------------------------------------------
  # Reuse candidate_observation for patients already processed.
  dbExecute(con, "
  CREATE OR REPLACE TABLE permissive_extra_ids AS
  SELECT p.patient_id, p.index_date
  FROM permissive_index0 p
  LEFT JOIN candidate_observation o USING (patient_id)
  WHERE o.patient_id IS NULL;
  ")

  extra_n <- dbGetQuery(
    con,
    "SELECT COUNT(*) AS n FROM permissive_extra_ids;"
  )$n[1]

  cat("\nNew permissive patients requiring encounter summarization:",
      extra_n, "\n")

  if (extra_n > 0) {

    if (
      REUSE_EXISTING_HEAVY_TABLES &&
      "permissive_extra_observation" %in% dbListTables(con)
    ) {
      cat("Reusing permissive_extra_observation.\n")
    } else {

      cat("\n============================================================\n")
      cat("Scanning encounter.csv for permissive extra patients...\n")
      cat("============================================================\n")
      t0 <- Sys.time()

      enc_sql <- paste0(
  "
  CREATE OR REPLACE TABLE permissive_extra_observation AS
  WITH enc AS (
    SELECT
      patient_id,
      TRY_CAST(STRPTIME(start_date, '%Y%m%d') AS DATE) AS encounter_date,
      type AS encounter_type
    FROM read_csv_auto(
      '", sql_path(ENCOUNTER_FILE), "',
      header=true,
      all_varchar=true,
      ignore_errors=false
    )
    WHERE start_date IS NOT NULL
      AND start_date <> ''
  )
  SELECT
    p.patient_id,

    MIN(e.encounter_date) AS first_encounter_date,
    MAX(e.encounter_date) AS last_encounter_date,

    COUNT(*) FILTER (
      WHERE e.encounter_date >= p.index_date - INTERVAL 365 DAY
        AND e.encounter_date < p.index_date
    ) AS n_encounters_prior_365,

    COUNT(*) FILTER (
      WHERE e.encounter_date >= p.index_date - INTERVAL 1095 DAY
        AND e.encounter_date < p.index_date
    ) AS n_encounters_prior_1095,

    COUNT(*) FILTER (
      WHERE e.encounter_date >= p.index_date - INTERVAL 365 DAY
        AND e.encounter_date < p.index_date
        AND e.encounter_type='AMB'
    ) AS n_ambulatory_prior_365,

    COUNT(*) FILTER (
      WHERE e.encounter_date >= p.index_date - INTERVAL 365 DAY
        AND e.encounter_date < p.index_date
        AND e.encounter_type='EMER'
    ) AS n_ed_prior_365,

    COUNT(*) FILTER (
      WHERE e.encounter_date >= p.index_date - INTERVAL 365 DAY
        AND e.encounter_date < p.index_date
        AND e.encounter_type IN ('IMP','ACUTE','NONAC','SS')
    ) AS n_inpatient_prior_365

  FROM permissive_extra_ids p
  LEFT JOIN enc e
    ON e.patient_id=p.patient_id

  GROUP BY p.patient_id;
  "
      )

      dbExecute(con, enc_sql)

      cat(
        "Extra encounter scan completed in ",
        round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1),
        " minutes.\n",
        sep=""
      )
    }
  }

  dbExecute(con, "
  CREATE OR REPLACE TABLE permissive_observation AS
  SELECT
    p.patient_id,
    o.first_encounter_date,
    o.last_encounter_date,
    o.n_encounters_prior_365,
    o.n_encounters_prior_1095,
    o.n_ambulatory_prior_365,
    o.n_ed_prior_365,
    o.n_inpatient_prior_365
  FROM permissive_index0 p
  JOIN candidate_observation o USING (patient_id)

  UNION ALL

  SELECT
    p.patient_id,
    o.first_encounter_date,
    o.last_encounter_date,
    o.n_encounters_prior_365,
    o.n_encounters_prior_1095,
    o.n_ambulatory_prior_365,
    o.n_ed_prior_365,
    o.n_inpatient_prior_365
  FROM permissive_index0 p
  JOIN permissive_extra_observation o USING (patient_id);
  ")

  # ------------------------------------------------------------
  # 4. DIAGNOSIS SCAN
  # ------------------------------------------------------------
  if (
    REUSE_EXISTING_HEAVY_TABLES &&
    "permissive_dx_summary" %in% dbListTables(con)
  ) {
    cat("\nReusing permissive_dx_summary; diagnosis.csv will NOT be rescanned.\n")
  } else {

    cat("\n============================================================\n")
    cat("Scanning diagnosis.csv once for permissive cohort...\n")
    cat("============================================================\n")
    t0 <- Sys.time()

    dx_template <- "
  CREATE OR REPLACE TABLE permissive_dx_summary AS

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
      AND UPPER(TRIM(code_system)) IN ('ICD-10-CM','ICD-9-CM')
  ),

  dx0 AS (
    SELECT *
    FROM dxraw
    WHERE
      (
        code_system='ICD-10-CM'
        AND (
             code_norm LIKE 'N40%'
          OR code_norm='C61'
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
          OR code_norm='N390'
          OR code_norm LIKE 'R31%'
          OR code_norm LIKE 'N21%'
        )
      )
      OR
      (
        code_system='ICD-9-CM'
        AND (
             code_norm LIKE '600%'
          OR code_norm='185'
          OR code_norm IN ('2252','2376','1921')
          OR code_norm LIKE '401%'
          OR code_norm LIKE '402%'
          OR code_norm LIKE '403%'
          OR code_norm LIKE '404%'
          OR code_norm LIKE '405%'
          OR code_norm LIKE '250%'
          OR code_norm LIKE '2780%'
          OR code_norm='278'
          OR code_norm LIKE '585%'
          OR code_norm LIKE '3051%'
          OR code_norm='V1582'
          OR code_norm LIKE '414%'
          OR code_norm LIKE '428%'
          OR (code_norm >= '430' AND code_norm < '439')
          OR code_norm LIKE '272%'
          OR code_norm LIKE '7330%'
          OR code_norm LIKE '704%'
          OR code_norm LIKE '7882%'
          OR code_norm='5990'
          OR code_norm LIKE '5997%'
          OR code_norm LIKE '594%'
        )
      )
  ),

  dx AS (
    SELECT
      d.*,

      CASE
        WHEN code_system='ICD-10-CM' AND code_norm LIKE 'N40%' THEN 1
        WHEN code_system='ICD-9-CM'  AND code_norm LIKE '600%' THEN 1
        ELSE 0
      END AS is_bph,

      CASE
        WHEN code_system='ICD-10-CM' AND code_norm='C61' THEN 1
        WHEN code_system='ICD-9-CM'  AND code_norm='185' THEN 1
        ELSE 0
      END AS is_prostate_cancer,

      CASE
        WHEN code_system='ICD-10-CM'
         AND code_norm IN ('D320','D329') THEN 1
        WHEN code_system='ICD-9-CM'
         AND code_norm='2252' THEN 1
        ELSE 0
      END AS is_strict_meningioma,

      CASE
        WHEN code_system='ICD-10-CM'
         AND code_norm='D320' THEN 1
        WHEN code_system='ICD-9-CM'
         AND code_norm='2252' THEN 1
        ELSE 0
      END AS is_specific_intracranial,

      CASE
        WHEN code_system='ICD-10-CM'
         AND code_norm='D320' THEN 1
        ELSE 0
      END AS is_d320,

      CASE
        WHEN code_system='ICD-10-CM'
         AND code_norm IN ('D320','D329','D420','D429','C700','C709')
          THEN 1
        WHEN code_system='ICD-9-CM'
         AND code_norm IN ('2252','2376','1921')
          THEN 1
        ELSE 0
      END AS is_broad_meningeal,

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
      p.patient_id,
      p.index_date,
      d.*
    FROM permissive_index0 p
    LEFT JOIN dx d
      ON d.patient_id=p.patient_id
  )

  SELECT
    patient_id,

    MAX(CASE
      WHEN is_bph=1
       AND dx_date >= index_date - INTERVAL {{BPH_LOOKBACK_DAYS}} DAY
       AND dx_date <= index_date + INTERVAL {{BPH_GRACE_DAYS}} DAY
      THEN 1 ELSE 0 END
    ) AS bph_near_index,

    MAX(CASE
      WHEN is_prostate_cancer=1
       AND dx_date < index_date
      THEN 1 ELSE 0 END
    ) AS prior_prostate_cancer,

    MIN(dx_date) FILTER (
      WHERE is_broad_meningeal=1
        AND dx_date <= index_date
    ) AS broad_preindex_date,

    MIN(dx_date) FILTER (
      WHERE is_broad_meningeal=1
        AND dx_date >= index_date
        AND dx_date < index_date + INTERVAL 365 DAY
    ) AS broad_lag1_date,

    MIN(dx_date) FILTER (
      WHERE is_strict_meningioma=1
        AND dx_date > index_date
    ) AS strict_event_post_index,

    MIN(dx_date) FILTER (
      WHERE is_strict_meningioma=1
        AND dx_date >= index_date + INTERVAL 365 DAY
    ) AS strict_event_post_lag1,

    MIN(dx_date) FILTER (
      WHERE is_specific_intracranial=1
        AND dx_date >= index_date + INTERVAL 365 DAY
    ) AS specific_event_post_lag1,

    MIN(dx_date) FILTER (
      WHERE is_d320=1
        AND dx_date >= index_date + INTERVAL 365 DAY
    ) AS d320_event_post_lag1,

    MAX(CASE WHEN is_hypertension=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS hypertension,

    MAX(CASE WHEN is_diabetes=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS diabetes,

    MAX(CASE WHEN is_obesity=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS obesity,

    MAX(CASE WHEN is_ckd=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS ckd,

    MAX(CASE WHEN is_smoking=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS smoking,

    MAX(CASE WHEN is_cad=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS cad,

    MAX(CASE WHEN is_hf=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS heart_failure,

    MAX(CASE WHEN is_cerebrovascular=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS cerebrovascular_disease,

    MAX(CASE WHEN is_hyperlipidemia=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS hyperlipidemia,

    MAX(CASE WHEN is_osteoporosis=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS osteoporosis,

    MAX(CASE WHEN is_alopecia=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS alopecia,

    MAX(CASE WHEN is_urinary_retention=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS urinary_retention,

    MAX(CASE WHEN is_uti=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS uti,

    MAX(CASE WHEN is_hematuria=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS hematuria,

    MAX(CASE WHEN is_bladder_stone=1
              AND dx_date >= index_date - INTERVAL {{BASELINE_DAYS}} DAY
              AND dx_date < index_date THEN 1 ELSE 0 END) AS bladder_stone

  FROM j
  GROUP BY patient_id;
  "

    dx_sql <- sql_fill(
      dx_template,
      list(
        DIAGNOSIS_FILE = sql_path(DIAGNOSIS_FILE),
        BPH_LOOKBACK_DAYS = BPH_LOOKBACK_DAYS,
        BPH_GRACE_DAYS = BPH_GRACE_DAYS,
        BASELINE_DAYS = BASELINE_DAYS
      )
    )

    dbExecute(con, dx_sql)

    cat(
      "Permissive diagnosis scan completed in ",
      round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1),
      " minutes.\n",
      sep=""
    )
  }

  # ------------------------------------------------------------
  # 5. MERGED PERMISSIVE DATASET
  # ------------------------------------------------------------
  dbExecute(con, "
  CREATE OR REPLACE TABLE permissive_eligibility AS
  SELECT
    p.*,
    s.index_finasteride_strength,
    s.index_dutasteride_strength,

    CASE
      WHEN p.first_day_5ari=1
       AND UPPER(TRIM(COALESCE(s.index_finasteride_strength,'')))='1 MG'
      THEN 1 ELSE 0
    END AS known_finasteride_1mg,

    o.first_encounter_date,
    o.last_encounter_date,
    o.n_encounters_prior_365,
    o.n_encounters_prior_1095,
    o.n_ambulatory_prior_365,
    o.n_ed_prior_365,
    o.n_inpatient_prior_365,

    COALESCE(d.bph_near_index,0) AS bph_near_index,
    COALESCE(d.prior_prostate_cancer,0) AS prior_prostate_cancer,
    d.broad_preindex_date,
    d.broad_lag1_date,
    d.strict_event_post_index,
    d.strict_event_post_lag1,
    d.specific_event_post_lag1,
    d.d320_event_post_lag1,

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
      WHEN p.exposure='5ARI'
       AND p.first_uro_alpha IS NOT NULL
       AND p.first_uro_alpha > p.index_date
       AND p.first_uro_alpha <= p.index_date + INTERVAL 30 DAY
        THEN 1
      WHEN p.exposure='URO_ALPHA'
       AND p.first_5ari IS NOT NULL
       AND p.first_5ari > p.index_date
       AND p.first_5ari <= p.index_date + INTERVAL 30 DAY
        THEN 1
      ELSE 0
    END AS opposite_class_within_30d

  FROM permissive_index0 p
  LEFT JOIN permissive_index_strength s USING (patient_id)
  LEFT JOIN permissive_observation o USING (patient_id)
  LEFT JOIN permissive_dx_summary d USING (patient_id);
  ")

  # ------------------------------------------------------------
  # 6. MAXIMALLY PERMISSIVE, 1-Y LAG
  # ------------------------------------------------------------
  dbExecute(con, paste0(
  "
  CREATE OR REPLACE TABLE analysis_cohort_permissive_max AS
  SELECT
    p.*,

    p.index_date + INTERVAL 365 DAY AS landmark_date,

    CASE
      WHEN p.strict_event_post_lag1 IS NOT NULL
       AND p.strict_event_post_lag1 <= p.last_encounter_date
      THEN 1 ELSE 0
    END AS meningioma_event_strict,

    CASE
      WHEN p.specific_event_post_lag1 IS NOT NULL
       AND p.specific_event_post_lag1 <= p.last_encounter_date
      THEN 1 ELSE 0
    END AS meningioma_event_specific,

    CASE
      WHEN p.strict_event_post_lag1 IS NOT NULL
       AND p.strict_event_post_lag1 <= p.last_encounter_date
      THEN p.strict_event_post_lag1
      ELSE p.last_encounter_date
    END AS followup_end_strict,

    DATE_DIFF(
      'day',
      p.index_date + INTERVAL 365 DAY,
      CASE
        WHEN p.strict_event_post_lag1 IS NOT NULL
         AND p.strict_event_post_lag1 <= p.last_encounter_date
        THEN p.strict_event_post_lag1
        ELSE p.last_encounter_date
      END
    ) AS followup_days_strict

  FROM permissive_eligibility p

  WHERE p.sex='M'
    AND p.age_at_index BETWEEN ", MIN_AGE, " AND ", MAX_AGE, "
    AND p.index_date >= DATE '", format(INDEX_START), "'
    AND p.index_date <= DATE '", format(INDEX_END_LAG1), "'
    AND p.first_encounter_date <= p.index_date - INTERVAL 365 DAY
    AND p.n_encounters_prior_365 >= 1
    AND p.broad_preindex_date IS NULL
    AND p.broad_lag1_date IS NULL
    AND p.last_encounter_date > p.index_date + INTERVAL 365 DAY;
  "
  ))

  # ------------------------------------------------------------
  # 7. PERMISSIVE BPH, 1-Y LAG
  # ------------------------------------------------------------
  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_permissive_bph AS
  SELECT *
  FROM analysis_cohort_permissive_max
  WHERE bph_near_index=1
    AND prior_prostate_cancer=0;
  ")

  # ------------------------------------------------------------
  # 8. MAXIMALLY PERMISSIVE ZERO-LAG
  # ------------------------------------------------------------
  dbExecute(con, paste0(
  "
  CREATE OR REPLACE TABLE analysis_cohort_permissive_zero_lag AS
  SELECT
    p.*,

    CASE
      WHEN p.strict_event_post_index IS NOT NULL
       AND p.strict_event_post_index <= p.last_encounter_date
      THEN 1 ELSE 0
    END AS meningioma_event_strict,

    CASE
      WHEN p.strict_event_post_index IS NOT NULL
       AND p.strict_event_post_index <= p.last_encounter_date
      THEN p.strict_event_post_index
      ELSE p.last_encounter_date
    END AS followup_end_strict,

    DATE_DIFF(
      'day',
      p.index_date,
      CASE
        WHEN p.strict_event_post_index IS NOT NULL
         AND p.strict_event_post_index <= p.last_encounter_date
        THEN p.strict_event_post_index
        ELSE p.last_encounter_date
      END
    ) AS followup_days_strict

  FROM permissive_eligibility p

  WHERE p.sex='M'
    AND p.age_at_index BETWEEN ", MIN_AGE, " AND ", MAX_AGE, "
    AND p.index_date >= DATE '", format(INDEX_START), "'
    AND p.index_date <= DATE '", format(INDEX_END_ZERO), "'
    AND p.first_encounter_date <= p.index_date - INTERVAL 365 DAY
    AND p.n_encounters_prior_365 >= 1
    AND p.broad_preindex_date IS NULL
    AND p.last_encounter_date > p.index_date + INTERVAL 30 DAY;
  "
  ))

  # ------------------------------------------------------------
  # 9. CLEAN ZERO-LAG LATENCY COHORT
  # ------------------------------------------------------------
  # This recovers the original clean active-comparator logic WITHOUT excluding
  # first-year outcomes and WITHOUT requiring one full year of future follow-up.
  # It is the preferred cohort for the interval-specific 0-1 / 1-2 / 2-5 / 5+
  # latency sensitivity.
  dbExecute(con, paste0(
  "
  CREATE OR REPLACE TABLE analysis_cohort_latency_clean AS
  SELECT
    p.*,

    CASE
      WHEN p.strict_event_post_index IS NOT NULL
       AND p.strict_event_post_index <= p.last_encounter_date
      THEN 1 ELSE 0
    END AS meningioma_event_strict,

    CASE
      WHEN p.strict_event_post_index IS NOT NULL
       AND p.strict_event_post_index <= p.last_encounter_date
      THEN p.strict_event_post_index
      ELSE p.last_encounter_date
    END AS followup_end_strict,

    DATE_DIFF(
      'day',
      p.index_date,
      CASE
        WHEN p.strict_event_post_index IS NOT NULL
         AND p.strict_event_post_index <= p.last_encounter_date
        THEN p.strict_event_post_index
        ELSE p.last_encounter_date
      END
    ) AS followup_days_strict

  FROM permissive_eligibility p

  WHERE p.sex='M'
    AND p.age_at_index BETWEEN ", MIN_AGE, " AND ", MAX_AGE, "
    AND p.index_date >= DATE '", format(INDEX_START), "'
    AND p.index_date <= DATE '", format(INDEX_END_ZERO), "'
    AND p.first_encounter_date <= p.index_date - INTERVAL 365 DAY
    AND p.n_encounters_prior_365 >= 1
    AND p.bph_near_index=1
    AND p.prior_prostate_cancer=0
    AND p.broad_preindex_date IS NULL
    AND p.n_drugs_on_first_day=1
    AND p.opposite_class_within_30d=0
    AND p.known_finasteride_1mg=0
    AND p.last_encounter_date > p.index_date + INTERVAL 30 DAY;
  "
  ))

  # ------------------------------------------------------------
  # 10. QC
  # ------------------------------------------------------------
  qc <- as.data.table(dbGetQuery(con, "
  SELECT 'permissive_max_1y_lag' AS cohort, exposure,
         COUNT(*) AS patients,
         SUM(meningioma_event_strict) AS events
  FROM analysis_cohort_permissive_max
  GROUP BY exposure

  UNION ALL

  SELECT 'permissive_bph_1y_lag' AS cohort, exposure,
         COUNT(*) AS patients,
         SUM(meningioma_event_strict) AS events
  FROM analysis_cohort_permissive_bph
  GROUP BY exposure

  UNION ALL

  SELECT 'permissive_zero_lag' AS cohort, exposure,
         COUNT(*) AS patients,
         SUM(meningioma_event_strict) AS events
  FROM analysis_cohort_permissive_zero_lag
  GROUP BY exposure

  UNION ALL

  SELECT 'latency_clean_zero_lag' AS cohort, exposure,
         COUNT(*) AS patients,
         SUM(meningioma_event_strict) AS events
  FROM analysis_cohort_latency_clean
  GROUP BY exposure

  ORDER BY cohort, exposure;
  "))

  first_day_qc <- as.data.table(dbGetQuery(con, "
  SELECT
    exposure,
    index_drug,
    COUNT(*) AS patients
  FROM analysis_cohort_permissive_max
  GROUP BY exposure, index_drug
  ORDER BY exposure, patients DESC;
  "))

  fwrite(qc, file.path(QC_DIR, "08_permissive_cohort_counts.csv"))
  fwrite(
    first_day_qc,
    file.path(QC_DIR, "08_permissive_first_day_regimens.csv")
  )

  cat("\n================ PERMISSIVE COHORT COUNTS ===========\n")
  print(qc)

  cat("\n================ TOP FIRST-DAY REGIMENS =============\n")
  print(first_day_qc[1:min(.N, 30)])

  cat("\nStage 08 complete.\n")

  dbDisconnect(con, shutdown = TRUE)
}

run_stage_08()
