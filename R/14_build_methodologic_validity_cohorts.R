# ============================================================
# 14_build_methodologic_validity_cohorts.R
#
# Finasteride / Dutasteride and Incident Meningioma
#
# Targeted methodological strengthening based on contemporary
# TriNetX / pharmacoepidemiology guidance.
#
# Builds:
#   1) Target-trial-aligned day-0 treatment cohorts
#      a. zero-lag ITT-like initiation effect
#      b. 1-y biologic-latency landmark
#   2) Prespecified negative-control outcome cohorts
#      a. benign cranial nerve neoplasm: D33.3 / ICD-9 225.1
#      b. benign spinal cord neoplasm: D33.4 / ICD-9 225.3
#   3) 180-day post-landmark utilization bins for IPCW
#   4) Internal meningioma phenotype-corroboration tables
#   5) TARGET-style target-trial specification table
#
# Heavy scans:
#   diagnosis.csv ONCE (checkpointed)
#   encounter.csv ONCE (checkpointed)
#
# IMPORTANT DESIGN DISTINCTION:
# The target-trial-aligned cohort determines initial treatment
# strictly from information available on day 0. It:
#   - requires BPH documented during the 365 d BEFORE/on index;
#   - excludes same-day 5ARI + uro-alpha cross-class starts;
#   - DOES NOT use treatment during days 1-30 to determine assignment;
#   - DOES NOT exclude later switching in the ITT-like analysis.
#
# The 1-y landmark is retained only as a biologic-latency analysis.
# ============================================================

run_stage_14 <- function() {

  options(stringsAsFactors = FALSE, scipen = 999)

  required_pkgs <- c("DBI", "duckdb", "data.table")
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0L) install.packages(missing_pkgs)

  suppressPackageStartupMessages({
    library(DBI)
    library(duckdb)
    library(data.table)
  })

  # ----------------------------------------------------------
  # 0. PATHS / SETTINGS
  # ----------------------------------------------------------
  RAW_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/Data/BPH"
  WORK_DIR <- "C:/TriNetX/bph_m_work"
  DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")
  TEMP_DIR <- file.path(WORK_DIR, "duckdb_temp")
  QC_DIR <- file.path(WORK_DIR, "qc_methodologic_validity")

  DIAGNOSIS_FILE <- file.path(RAW_DIR, "diagnosis.csv")
  ENCOUNTER_FILE <- file.path(RAW_DIR, "encounter.csv")

  dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

  N_THREADS <- min(
    24L,
    max(4L, parallel::detectCores(logical = TRUE) - 6L)
  )
  MEMORY_LIMIT <- "48GB"
  REUSE_EXISTING_HEAVY_TABLES <- TRUE

  MIN_AGE <- 40L
  MAX_AGE <- 89L
  INDEX_START <- as.Date("2004-01-01")
  BPH_PREINDEX_DAYS <- 365L
  IPCW_INTERVAL_DAYS <- 180L

  sql_path <- function(x) {
    x <- normalizePath(x, winslash = "/", mustWork = FALSE)
    gsub("'", "''", x, fixed = TRUE)
  }

  sql_date <- function(x) format(as.Date(x), "%Y-%m-%d")

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
  dbExecute(
    con,
    sprintf(
      "PRAGMA temp_directory='%s';",
      sql_path(TEMP_DIR)
    )
  )
  dbExecute(con, "SET preserve_insertion_order=false;")

  required_tables <- c(
    "permissive_index0",
    "permissive_eligibility",
    "candidate_observation",
    "analysis_cohort_primary_revised"
  )

  missing_tables <- setdiff(
    required_tables,
    dbListTables(con)
  )

  if (length(missing_tables) > 0L) {
    stop(
      "Missing prerequisite tables. Run Stages 05-12 first. Missing: ",
      paste(missing_tables, collapse = ", ")
    )
  }

  data_cutoff <- as.Date(
    dbGetQuery(
      con,
      "SELECT MAX(last_encounter_date) AS x
       FROM candidate_observation;"
    )$x[1]
  )

  if (is.na(data_cutoff)) {
    stop("Could not derive observed data cutoff.")
  }

  index_end_zero <- data_cutoff - 1L
  index_end_lag1 <- data_cutoff - 365L

  cat("\nObserved data cutoff:", format(data_cutoff), "\n")
  cat("Day-0 ITT-like index end:", format(index_end_zero), "\n")
  cat("1-y landmark index end:", format(index_end_lag1), "\n")

  # ----------------------------------------------------------
  # 1. ONE DIAGNOSIS SCAN
  # ----------------------------------------------------------
  #
  # In addition to the meningioma phenotype, this scan obtains:
  #   - exact pre-index BPH timing (no +30 d grace);
  #   - prostate cancer;
  #   - two prespecified neurologic negative controls.
  #
  # Negative controls:
  #   D33.3 / 225.1 = benign neoplasm of cranial nerves
  #   D33.4 / 225.3 = benign neoplasm of spinal cord
  #
  if (
    REUSE_EXISTING_HEAVY_TABLES &&
    "method_dx_events" %in% dbListTables(con)
  ) {

    cat(
      "\nReusing method_dx_events; diagnosis.csv will NOT be rescanned.\n"
    )

  } else {

    cat("\n============================================================\n")
    cat("Scanning diagnosis.csv for methodological-validity phenotypes...\n")
    cat("============================================================\n")
    t0 <- Sys.time()

    dx_sql <- paste0(
"
CREATE OR REPLACE TABLE method_dx_events AS

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
    AND UPPER(TRIM(code_system)) IN ('ICD-10-CM','ICD-9-CM')
),

f AS (
  SELECT
    patient_id,
    code_system,
    code_norm,
    dx_date,

    CASE
      WHEN code_system='ICD-10-CM'
       AND code_norm LIKE 'N40%' THEN 1
      WHEN code_system='ICD-9-CM'
       AND code_norm LIKE '600%' THEN 1
      ELSE 0
    END AS is_bph,

    CASE
      WHEN code_system='ICD-10-CM'
       AND code_norm='C61' THEN 1
      WHEN code_system='ICD-9-CM'
       AND code_norm='185' THEN 1
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
       AND code_norm IN ('D320','D329','D420','D429','C700','C709')
        THEN 1
      WHEN code_system='ICD-9-CM'
       AND code_norm IN ('2252','2376','1921')
        THEN 1
      ELSE 0
    END AS is_broad_meningeal,

    CASE
      WHEN code_system='ICD-10-CM'
       AND code_norm='D333' THEN 1
      WHEN code_system='ICD-9-CM'
       AND code_norm='2251' THEN 1
      ELSE 0
    END AS is_nco_cranial_nerve,

    CASE
      WHEN code_system='ICD-10-CM'
       AND code_norm='D334' THEN 1
      WHEN code_system='ICD-9-CM'
       AND code_norm='2253' THEN 1
      ELSE 0
    END AS is_nco_spinal_cord

  FROM d

  WHERE
    (
      code_system='ICD-10-CM'
      AND (
           code_norm LIKE 'N40%'
        OR code_norm='C61'
        OR code_norm IN (
          'D320','D329','D420','D429','C700','C709',
          'D333','D334'
        )
      )
    )
    OR
    (
      code_system='ICD-9-CM'
      AND (
           code_norm LIKE '600%'
        OR code_norm='185'
        OR code_norm IN (
          '2252','2376','1921','2251','2253'
        )
      )
    )
)

SELECT DISTINCT
  f.*
FROM f
INNER JOIN permissive_index0 p
  ON p.patient_id=f.patient_id
WHERE f.dx_date IS NOT NULL;
"
    )

    dbExecute(con, dx_sql)

    try(
      dbExecute(
        con,
        "CREATE INDEX IF NOT EXISTS idx_method_dx_patient_date
         ON method_dx_events(patient_id, dx_date);"
      ),
      silent = TRUE
    )

    cat(
      "Diagnosis scan completed in ",
      round(
        as.numeric(
          difftime(Sys.time(), t0, units = "mins")
        ),
        1
      ),
      " minutes.\n",
      sep = ""
    )
  }

  # ----------------------------------------------------------
  # 2. TEMPORAL DIAGNOSIS SUMMARY
  # ----------------------------------------------------------
  dbExecute(
    con,
"
CREATE OR REPLACE TABLE method_dx_summary AS

WITH j AS (
  SELECT
    p.patient_id,
    p.index_date,
    d.*
  FROM permissive_index0 p
  LEFT JOIN method_dx_events d
    ON d.patient_id=p.patient_id
)

SELECT
  patient_id,

  MAX(
    CASE
      WHEN is_bph=1
       AND dx_date >= index_date - INTERVAL 365 DAY
       AND dx_date <= index_date
      THEN 1 ELSE 0
    END
  ) AS bph_preindex_365,

  MAX(
    CASE
      WHEN is_prostate_cancer=1
       AND dx_date < index_date
      THEN 1 ELSE 0
    END
  ) AS prior_prostate_cancer,

  MIN(dx_date) FILTER (
    WHERE is_broad_meningeal=1
      AND dx_date <= index_date
  ) AS broad_preindex_or_same_day,

  MIN(dx_date) FILTER (
    WHERE is_broad_meningeal=1
      AND dx_date > index_date
      AND dx_date < index_date + INTERVAL 365 DAY
  ) AS broad_lag1_date,

  MIN(dx_date) FILTER (
    WHERE is_strict_meningioma=1
      AND dx_date > index_date
  ) AS strict_post_index_first,

  MIN(dx_date) FILTER (
    WHERE is_strict_meningioma=1
      AND dx_date >= index_date + INTERVAL 365 DAY
  ) AS strict_post_lag1_first,

  MIN(dx_date) FILTER (
    WHERE is_specific_intracranial=1
      AND dx_date >= index_date + INTERVAL 365 DAY
  ) AS specific_post_lag1_first,

  MIN(dx_date) FILTER (
    WHERE is_nco_cranial_nerve=1
      AND dx_date < index_date + INTERVAL 365 DAY
  ) AS nco_cranial_pre_landmark,

  MIN(dx_date) FILTER (
    WHERE is_nco_cranial_nerve=1
      AND dx_date >= index_date + INTERVAL 365 DAY
  ) AS nco_cranial_post_landmark,

  MIN(dx_date) FILTER (
    WHERE is_nco_spinal_cord=1
      AND dx_date < index_date + INTERVAL 365 DAY
  ) AS nco_spinal_pre_landmark,

  MIN(dx_date) FILTER (
    WHERE is_nco_spinal_cord=1
      AND dx_date >= index_date + INTERVAL 365 DAY
  ) AS nco_spinal_post_landmark

FROM j
GROUP BY patient_id;
"
  )

  # ----------------------------------------------------------
  # 3. TARGET-TRIAL-ALIGNED DAY-0 COHORT: ZERO LAG
  # ----------------------------------------------------------
  #
  # Treatment is determined by day-0 class only.
  # Later switching is ignored (ITT-like initiation strategy).
  #
  dbExecute(
    con,
    paste0(
"
CREATE OR REPLACE TABLE analysis_cohort_tte_index_itt_zero AS

SELECT
  p.*,

  m.bph_preindex_365,
  m.broad_preindex_or_same_day,
  m.strict_post_index_first,

  CASE
    WHEN m.strict_post_index_first IS NOT NULL
     AND m.strict_post_index_first <= p.last_encounter_date
    THEN 1 ELSE 0
  END AS meningioma_event_strict,

  CASE
    WHEN m.strict_post_index_first IS NOT NULL
     AND m.strict_post_index_first <= p.last_encounter_date
    THEN m.strict_post_index_first
    ELSE p.last_encounter_date
  END AS followup_end_strict,

  DATE_DIFF(
    'day',
    p.index_date,
    CASE
      WHEN m.strict_post_index_first IS NOT NULL
       AND m.strict_post_index_first <= p.last_encounter_date
      THEN m.strict_post_index_first
      ELSE p.last_encounter_date
    END
  ) AS followup_days_strict

FROM permissive_eligibility p
JOIN method_dx_summary m USING (patient_id)

WHERE p.exposure IN ('5ARI','URO_ALPHA')
  AND NOT (
    p.first_day_5ari=1
    AND p.first_day_uro_alpha=1
  )
  AND p.sex='M'
  AND p.age_at_index BETWEEN ", MIN_AGE, " AND ", MAX_AGE, "
  AND p.index_date >= DATE '", sql_date(INDEX_START), "'
  AND p.index_date <= DATE '", sql_date(index_end_zero), "'
  AND p.first_encounter_date <= p.index_date - INTERVAL 365 DAY
  AND p.n_encounters_prior_365 >= 1
  AND COALESCE(m.bph_preindex_365,0)=1
  AND COALESCE(m.prior_prostate_cancer,0)=0
  AND m.broad_preindex_or_same_day IS NULL
  AND COALESCE(p.known_finasteride_1mg,0)=0
  AND p.last_encounter_date > p.index_date;
"
    )
  )

  # ----------------------------------------------------------
  # 4. TARGET-TRIAL-ALIGNED DAY-0 COHORT: 1-Y LANDMARK
  # ----------------------------------------------------------
  dbExecute(
    con,
    paste0(
"
CREATE OR REPLACE TABLE analysis_cohort_tte_index_itt_lag1 AS

SELECT
  p.*,

  m.bph_preindex_365,
  m.broad_preindex_or_same_day,
  m.broad_lag1_date,
  m.strict_post_lag1_first,

  p.index_date + INTERVAL 365 DAY AS landmark_date,

  CASE
    WHEN m.strict_post_lag1_first IS NOT NULL
     AND m.strict_post_lag1_first <= p.last_encounter_date
    THEN 1 ELSE 0
  END AS meningioma_event_strict,

  CASE
    WHEN m.strict_post_lag1_first IS NOT NULL
     AND m.strict_post_lag1_first <= p.last_encounter_date
    THEN m.strict_post_lag1_first
    ELSE p.last_encounter_date
  END AS followup_end_strict,

  DATE_DIFF(
    'day',
    p.index_date + INTERVAL 365 DAY,
    CASE
      WHEN m.strict_post_lag1_first IS NOT NULL
       AND m.strict_post_lag1_first <= p.last_encounter_date
      THEN m.strict_post_lag1_first
      ELSE p.last_encounter_date
    END
  ) AS followup_days_strict

FROM permissive_eligibility p
JOIN method_dx_summary m USING (patient_id)

WHERE p.exposure IN ('5ARI','URO_ALPHA')
  AND NOT (
    p.first_day_5ari=1
    AND p.first_day_uro_alpha=1
  )
  AND p.sex='M'
  AND p.age_at_index BETWEEN ", MIN_AGE, " AND ", MAX_AGE, "
  AND p.index_date >= DATE '", sql_date(INDEX_START), "'
  AND p.index_date <= DATE '", sql_date(index_end_lag1), "'
  AND p.first_encounter_date <= p.index_date - INTERVAL 365 DAY
  AND p.n_encounters_prior_365 >= 1
  AND COALESCE(m.bph_preindex_365,0)=1
  AND COALESCE(m.prior_prostate_cancer,0)=0
  AND m.broad_preindex_or_same_day IS NULL
  AND m.broad_lag1_date IS NULL
  AND COALESCE(p.known_finasteride_1mg,0)=0
  AND p.last_encounter_date > p.index_date + INTERVAL 365 DAY;
"
    )
  )

  # ----------------------------------------------------------
  # 5. NEGATIVE-CONTROL OUTCOME COHORTS
  # ----------------------------------------------------------
  #
  # These use the FINAL primary treatment cohort and its 1-y risk
  # landmark, but substitute outcomes that should not plausibly be
  # caused by 5-alpha-reductase inhibition.
  #
  dbExecute(
    con,
"
CREATE OR REPLACE TABLE analysis_cohort_nco_cranial_nerve AS

SELECT
  a.*,

  CASE
    WHEN m.nco_cranial_post_landmark IS NOT NULL
     AND m.nco_cranial_post_landmark <= a.last_encounter_date
    THEN 1 ELSE 0
  END AS nco_event,

  CASE
    WHEN m.nco_cranial_post_landmark IS NOT NULL
     AND m.nco_cranial_post_landmark <= a.last_encounter_date
    THEN m.nco_cranial_post_landmark
    ELSE a.last_encounter_date
  END AS nco_followup_end,

  DATE_DIFF(
    'day',
    a.landmark_date,
    CASE
      WHEN m.nco_cranial_post_landmark IS NOT NULL
       AND m.nco_cranial_post_landmark <= a.last_encounter_date
      THEN m.nco_cranial_post_landmark
      ELSE a.last_encounter_date
    END
  ) AS nco_followup_days

FROM analysis_cohort_primary_revised a
JOIN method_dx_summary m USING (patient_id)

WHERE m.nco_cranial_pre_landmark IS NULL
  AND a.last_encounter_date > a.landmark_date;
"
  )

  dbExecute(
    con,
"
CREATE OR REPLACE TABLE analysis_cohort_nco_spinal_cord AS

SELECT
  a.*,

  CASE
    WHEN m.nco_spinal_post_landmark IS NOT NULL
     AND m.nco_spinal_post_landmark <= a.last_encounter_date
    THEN 1 ELSE 0
  END AS nco_event,

  CASE
    WHEN m.nco_spinal_post_landmark IS NOT NULL
     AND m.nco_spinal_post_landmark <= a.last_encounter_date
    THEN m.nco_spinal_post_landmark
    ELSE a.last_encounter_date
  END AS nco_followup_end,

  DATE_DIFF(
    'day',
    a.landmark_date,
    CASE
      WHEN m.nco_spinal_post_landmark IS NOT NULL
       AND m.nco_spinal_post_landmark <= a.last_encounter_date
      THEN m.nco_spinal_post_landmark
      ELSE a.last_encounter_date
    END
  ) AS nco_followup_days

FROM analysis_cohort_primary_revised a
JOIN method_dx_summary m USING (patient_id)

WHERE m.nco_spinal_pre_landmark IS NULL
  AND a.last_encounter_date > a.landmark_date;
"
  )

  # ----------------------------------------------------------
  # 6. 180-DAY UTILIZATION BINS FOR IPCW
  # ----------------------------------------------------------
  #
  # Bin -1 = 180 d immediately BEFORE the 1-y risk landmark.
  # Bin  0 = first 180 d AFTER landmark, etc.
  #
  # The Stage-15 censoring model uses utilization from the
  # immediately PRECEDING bin, so all time-updated predictors are
  # known at the beginning of the weighted risk interval.
  #
  if (
    REUSE_EXISTING_HEAVY_TABLES &&
    "method_primary_180d_util" %in% dbListTables(con)
  ) {

    cat(
      "\nReusing method_primary_180d_util; encounter.csv will NOT be rescanned.\n"
    )

  } else {

    cat("\n============================================================\n")
    cat("Scanning encounter.csv for 180-d EHR-continuity utilization bins...\n")
    cat("============================================================\n")
    t0 <- Sys.time()

    enc_sql <- paste0(
"
CREATE OR REPLACE TABLE method_primary_180d_util AS

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
    a.patient_id,
    a.landmark_date,
    a.last_encounter_date,
    e.encounter_date,
    e.encounter_type,

    CAST(
      FLOOR(
        DATE_DIFF('day', a.landmark_date, e.encounter_date) /
        ", IPCW_INTERVAL_DAYS, ".0
      ) AS INTEGER
    ) AS interval_bin

  FROM analysis_cohort_primary_revised a
  JOIN e
    ON e.patient_id=a.patient_id
   AND e.encounter_date >=
       a.landmark_date - INTERVAL ", IPCW_INTERVAL_DAYS, " DAY
   AND e.encounter_date <= a.last_encounter_date
)

SELECT
  patient_id,
  interval_bin,

  COUNT(*) AS encounters_total,
  COUNT(DISTINCT encounter_date) AS encounter_days,

  COUNT(*) FILTER (
    WHERE encounter_type='AMB'
  ) AS ambulatory_encounters,

  COUNT(*) FILTER (
    WHERE encounter_type='EMER'
  ) AS ed_encounters,

  COUNT(*) FILTER (
    WHERE encounter_type IN ('IMP','ACUTE','NONAC','SS')
  ) AS inpatient_encounters

FROM j
GROUP BY patient_id, interval_bin;
"
    )

    dbExecute(con, enc_sql)

    try(
      dbExecute(
        con,
        "CREATE INDEX IF NOT EXISTS idx_method_primary_180d_util
         ON method_primary_180d_util(patient_id, interval_bin);"
      ),
      silent = TRUE
    )

    cat(
      "Encounter scan completed in ",
      round(
        as.numeric(
          difftime(Sys.time(), t0, units = "mins")
        ),
        1
      ),
      " minutes.\n",
      sep = ""
    )
  }

  # ----------------------------------------------------------
  # 7. PHENOTYPE CORROBORATION
  # ----------------------------------------------------------
  #
  # This is NOT claimed to be a formal PPV validation study.
  # It reports concordance of diagnosis-defined incident events
  # with independent / more-specific evidence available internally.
  #
  events <- as.data.table(
    dbGetQuery(
      con,
"
SELECT
  patient_id,
  exposure,
  index_date,
  landmark_date,
  last_encounter_date,
  incident_meningioma_date_strict,
  meningioma_confirmation_date
FROM analysis_cohort_primary_revised
WHERE meningioma_event_strict=1;
"
    )
  )

  events[
    ,
    `:=`(
      index_date = as.Date(index_date),
      landmark_date = as.Date(landmark_date),
      last_encounter_date = as.Date(last_encounter_date),
      incident_meningioma_date_strict =
        as.Date(incident_meningioma_date_strict),
      meningioma_confirmation_date =
        as.Date(meningioma_confirmation_date)
    )
  ]

  dx_first <- as.data.table(
    dbGetQuery(
      con,
"
SELECT
  patient_id,
  dx_date,
  MAX(is_specific_intracranial) AS high_specificity,
  MAX(
    CASE
      WHEN code_system='ICD-10-CM'
       AND code_norm='D329'
      THEN 1 ELSE 0
    END
  ) AS has_d329
FROM method_dx_events
WHERE is_strict_meningioma=1
GROUP BY patient_id, dx_date;
"
    )
  )

  dx_first[, dx_date := as.Date(dx_date)]

  events <- merge(
    events,
    dx_first,
    by.x = c(
      "patient_id",
      "incident_meningioma_date_strict"
    ),
    by.y = c(
      "patient_id",
      "dx_date"
    ),
    all.x = TRUE
  )

  events[
    ,
    `:=`(
      high_specificity = fifelse(
        is.na(high_specificity),
        0L,
        as.integer(high_specificity)
      ),
      first_dx_has_d329 = fifelse(
        is.na(has_d329),
        0L,
        as.integer(has_d329)
      ),
      repeat_code_ge30d = as.integer(
        !is.na(meningioma_confirmation_date) &
        meningioma_confirmation_date <= last_encounter_date
      ),
      tumor_table = 0L,
      cancer_registry = 0L,
      resection_or_srs = 0L,
      oncology_treatment_table = 0L
    )
  ]

  tbls <- dbListTables(con)

  # ICD-O tumor-table support within +/-365 d of the incident diagnosis.
  if ("tumor_meningioma_events" %in% tbls) {

    tumor <- as.data.table(
      dbGetQuery(
        con,
"
SELECT
  patient_id,
  diagnosis_date,
  source_id
FROM tumor_meningioma_events;
"
      )
    )

    tumor[, diagnosis_date := as.Date(diagnosis_date)]

    tmp <- merge(
      events[
        ,
        .(
          patient_id,
          event_date = incident_meningioma_date_strict
        )
      ],
      tumor,
      by = "patient_id",
      allow.cartesian = TRUE
    )

    tmp <- tmp[
      !is.na(diagnosis_date) &
      abs(
        as.numeric(diagnosis_date - event_date)
      ) <= 365
    ]

    if (nrow(tmp) > 0L) {
      tumor_ids <- unique(tmp$patient_id)
      registry_ids <- unique(
        tmp[
          grepl(
            "CANCER",
            toupper(fifelse(is.na(source_id), "", source_id))
          ),
          patient_id
        ]
      )

      events[
        patient_id %in% tumor_ids,
        tumor_table := 1L
      ]

      events[
        patient_id %in% registry_ids,
        cancer_registry := 1L
      ]
    }
  }

  # Meningioma-directed resection / SRS in the first 180 d after diagnosis.
  if ("meningioma_treatment_procedures" %in% tbls) {

    proc <- as.data.table(
      dbGetQuery(
        con,
"
SELECT
  patient_id,
  procedure_date,
  treatment_type
FROM meningioma_treatment_procedures;
"
      )
    )

    proc[, procedure_date := as.Date(procedure_date)]

    tmp <- merge(
      events[
        ,
        .(
          patient_id,
          event_date = incident_meningioma_date_strict
        )
      ],
      proc,
      by = "patient_id",
      allow.cartesian = TRUE
    )

    tmp <- tmp[
      !is.na(procedure_date) &
      procedure_date >= event_date &
      procedure_date <= event_date + 180L
    ]

    if (nrow(tmp) > 0L) {
      tx_ids <- unique(tmp$patient_id)

      events[
        patient_id %in% tx_ids,
        resection_or_srs := 1L
      ]
    }
  }

  # Oncology-treatment table support in the first year after diagnosis.
  if ("oncology_treated_meningioma_events" %in% tbls) {

    onc <- as.data.table(
      dbGetQuery(
        con,
"
SELECT
  patient_id,
  diagnosis_date,
  treatment_start_date
FROM oncology_treated_meningioma_events;
"
      )
    )

    onc[
      ,
      `:=`(
        diagnosis_date = as.Date(diagnosis_date),
        treatment_start_date = as.Date(treatment_start_date)
      )
    ]

    tmp <- merge(
      events[
        ,
        .(
          patient_id,
          event_date = incident_meningioma_date_strict
        )
      ],
      onc,
      by = "patient_id",
      allow.cartesian = TRUE
    )

    tmp <- tmp[
      !is.na(treatment_start_date) &
      treatment_start_date >= event_date &
      treatment_start_date <= event_date + 365L
    ]

    if (nrow(tmp) > 0L) {
      onc_ids <- unique(tmp$patient_id)

      events[
        patient_id %in% onc_ids,
        oncology_treatment_table := 1L
      ]
    }
  }

  events[
    ,
    independent_corroboration := as.integer(
      tumor_table == 1L |
      cancer_registry == 1L |
      resection_or_srs == 1L |
      oncology_treatment_table == 1L
    )
  ]

  events[
    ,
    any_corroboration := as.integer(
      high_specificity == 1L |
      repeat_code_ge30d == 1L |
      independent_corroboration == 1L
    )
  ]

  events[
    ,
    isolated_unspecified_single_code := as.integer(
      high_specificity == 0L &
      first_dx_has_d329 == 1L &
      repeat_code_ge30d == 0L &
      independent_corroboration == 0L
    )
  ]

  fwrite(
    events,
    file.path(
      QC_DIR,
      "14_phenotype_corroboration_patient_level.csv"
    )
  )

  criteria <- c(
    high_specificity = "First incident diagnosis included D32.0 / ICD-9 225.2",
    repeat_code_ge30d = "Second meningioma diagnosis >=30 d later",
    tumor_table = "ICD-O 953x tumor-table corroboration within +/-365 d",
    cancer_registry = "Cancer Registry corroboration within +/-365 d",
    resection_or_srs = "Meningioma resection / cranial SRS within 180 d",
    oncology_treatment_table = "Oncology-treatment-table corroboration within 365 d",
    independent_corroboration = "At least 1 independent tumor/treatment source",
    any_corroboration = "At least 1 corroborating criterion",
    isolated_unspecified_single_code = "Isolated D32.9 single-code event without corroboration"
  )

  phenotype_summary <- rbindlist(
    lapply(
      names(criteria),
      function(v) {

        overall_n <- sum(events[[v]] == 1L, na.rm = TRUE)
        overall_N <- nrow(events)

        e5 <- events[exposure == "5ARI"]
        eu <- events[exposure == "URO_ALPHA"]

        data.table(
          criterion = criteria[[v]],
          overall_n = overall_n,
          overall_N = overall_N,
          overall_pct = 100 * overall_n / overall_N,

          n_5ari = sum(e5[[v]] == 1L, na.rm = TRUE),
          N_5ari = nrow(e5),
          pct_5ari = 100 * mean(
            e5[[v]] == 1L,
            na.rm = TRUE
          ),

          n_uro_alpha = sum(eu[[v]] == 1L, na.rm = TRUE),
          N_uro_alpha = nrow(eu),
          pct_uro_alpha = 100 * mean(
            eu[[v]] == 1L,
            na.rm = TRUE
          )
        )
      }
    )
  )

  fwrite(
    phenotype_summary,
    file.path(
      QC_DIR,
      "14_phenotype_corroboration_summary.csv"
    )
  )

  # ----------------------------------------------------------
  # 8. TARGET-TRIAL SPECIFICATION TABLE
  # ----------------------------------------------------------
  target_trial <- data.table(
    component = c(
      "Eligibility",
      "Treatment strategies",
      "Treatment assignment",
      "Time zero",
      "Grace period",
      "Primary causal contrast",
      "Biologic-latency contrast",
      "Outcome",
      "Follow-up",
      "Confounding adjustment",
      "Competing event",
      "Loss to observable EHR follow-up",
      "Treatment switching"
    ),

    target_trial = c(
      "Men with BPH, no prior prostate cancer or meningeal neoplasm, adequate prior observation",
      "Initiate a 5-alpha-reductase inhibitor vs initiate a uroselective alpha-blocker",
      "Random assignment at treatment initiation",
      "Eligible treatment initiation",
      "None",
      "Assignment to initiation strategy regardless of subsequent switching",
      "Same initiation strategy with risk follow-up beginning at 1 year",
      "Incident meningioma",
      "Until outcome, death/loss of observation, administrative end, or specified horizon",
      "Randomization",
      "Death",
      "Administrative follow-up",
      "Ignored under intention-to-treat"
    ),

    observational_emulation = c(
      "Male, age 40-89 y; BPH documented in prior 365 d/on index; >=365 d observable history; no prior prostate cancer or broad meningeal tumor; known finasteride 1-mg initiators excluded",
      "First observed day-0 5ARI vs day-0 uroselective alpha-blocker; same-day cross-class starts excluded",
      "Observed treatment; propensity-score overlap weighting for measured confounding",
      "First eligible study-drug initiation",
      "None for treatment assignment; no days 1-30 treatment information used to define assignment",
      "ITT-like initiation effect beginning immediately after index",
      "Separate 1-y landmark sensitivity; lag-period meningeal tumors excluded before risk-set entry",
      "Diagnosis-defined D32.0/D32.9/ICD-9 225.2; high-specificity and registry/treatment definitions assessed separately",
      "Outcome or last observed encounter; administrative sensitivity analyses as specified",
      "Overlap weighting; balance assessed by standardized mean differences",
      "Fine-Gray sensitivity using month-year death marker",
      "IPCW sensitivity using 180-d person-period EHR continuity model",
      "Ignored in ITT-like analysis; sustained/persistent exposure analyses reported separately"
    )
  )

  fwrite(
    target_trial,
    file.path(
      QC_DIR,
      "14_target_trial_specification.csv"
    )
  )

  # ----------------------------------------------------------
  # 9. QC OUTPUTS
  # ----------------------------------------------------------
  tte_counts <- as.data.table(
    dbGetQuery(
      con,
"
SELECT
  'Day-0 ITT-like, zero lag' AS cohort,
  exposure,
  COUNT(*) AS patients,
  SUM(meningioma_event_strict) AS events,
  MEDIAN(followup_days_strict)/365.25 AS median_followup_y
FROM analysis_cohort_tte_index_itt_zero
GROUP BY exposure

UNION ALL

SELECT
  'Day-0 ITT-like, 1-y landmark' AS cohort,
  exposure,
  COUNT(*) AS patients,
  SUM(meningioma_event_strict) AS events,
  MEDIAN(followup_days_strict)/365.25 AS median_followup_y
FROM analysis_cohort_tte_index_itt_lag1
GROUP BY exposure

ORDER BY cohort, exposure;
"
    )
  )

  nco_counts <- as.data.table(
    dbGetQuery(
      con,
"
SELECT
  'Benign cranial nerve neoplasm (D33.3/225.1)' AS outcome,
  exposure,
  COUNT(*) AS patients,
  SUM(nco_event) AS events,
  MEDIAN(nco_followup_days)/365.25 AS median_followup_y
FROM analysis_cohort_nco_cranial_nerve
GROUP BY exposure

UNION ALL

SELECT
  'Benign spinal cord neoplasm (D33.4/225.3)' AS outcome,
  exposure,
  COUNT(*) AS patients,
  SUM(nco_event) AS events,
  MEDIAN(nco_followup_days)/365.25 AS median_followup_y
FROM analysis_cohort_nco_spinal_cord
GROUP BY exposure

ORDER BY outcome, exposure;
"
    )
  )

  util_qc <- as.data.table(
    dbGetQuery(
      con,
"
SELECT
  a.exposure,
  u.interval_bin,
  COUNT(DISTINCT u.patient_id) AS patients_with_activity,
  MEDIAN(u.encounter_days) AS median_encounter_days,
  MEDIAN(u.ambulatory_encounters) AS median_ambulatory,
  MEDIAN(u.ed_encounters) AS median_ed,
  MEDIAN(u.inpatient_encounters) AS median_inpatient
FROM method_primary_180d_util u
JOIN analysis_cohort_primary_revised a
  ON a.patient_id=u.patient_id
WHERE u.interval_bin IN (-1,0,1,2)
GROUP BY a.exposure, u.interval_bin
ORDER BY a.exposure, u.interval_bin;
"
    )
  )

  fwrite(
    tte_counts,
    file.path(QC_DIR, "14_tte_cohort_counts.csv")
  )

  fwrite(
    nco_counts,
    file.path(QC_DIR, "14_negative_control_counts.csv")
  )

  fwrite(
    util_qc,
    file.path(QC_DIR, "14_ipcw_utilization_qc.csv")
  )

  cat("\n================ TARGET-TRIAL COHORTS ===============\n")
  print(tte_counts)

  cat("\n================ NEGATIVE CONTROLS ===================\n")
  print(nco_counts)

  cat("\n================ PHENOTYPE CORROBORATION =============\n")
  print(phenotype_summary)

  cat("\nStage 14 complete.\n")
  cat("QC directory:", QC_DIR, "\n")
}

run_stage_14()
