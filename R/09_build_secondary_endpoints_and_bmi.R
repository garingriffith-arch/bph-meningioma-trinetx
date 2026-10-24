# NOTE: This stage is wrapped in a function so on.exit() is tied to
# the stage function rather than source()/eval(). This prevents the
# DuckDB connection from being closed immediately when sourced.

run_stage_09 <- function() {
  # ============================================================
  # 09_build_secondary_endpoints_and_bmi.R
  #
  # Builds independent / high-clinical-specificity endpoints:
  #   1) tumor-table ICD-O meningioma endpoint
  #   2) cancer-registry-only ICD-O endpoint
  #   3) meningioma diagnosis followed by meningioma-specific resection
  #      or cranial stereotactic radiosurgery within 180 days
  #   4) oncology-treatment-table meningioma endpoint, if available
  #   5) optional measured pre-index BMI from vitals_signs.csv
  #
  # New raw scans are checkpointed.
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

  TUMOR_FILE <- file.path(RAW_DIR, "tumor.csv")
  PROCEDURE_FILE <- file.path(RAW_DIR, "procedure.csv")
  ONCOLOGY_TREATMENT_FILE <- file.path(RAW_DIR, "oncology_treatment.csv")
  VITALS_FILE <- file.path(RAW_DIR, "vitals_signs.csv")

  # The BMI scan is ~78 GB in this export. It is enabled because you asked to
  # run the full literature-informed sensitivity package. Set FALSE if you want
  # to defer it without affecting the other endpoints.
  RUN_BMI_SCAN <- TRUE

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
    "meningioma_dx_events_revised",
    "standardized_terminology"
  )
  missing <- setdiff(needed, dbListTables(con))
  if (length(missing) > 0L) {
    stop("Missing required tables: ", paste(missing, collapse=", "))
  }

  # ------------------------------------------------------------
  # 1. TUMOR-TABLE MENINGIOMA
  # ------------------------------------------------------------
  if (file.exists(TUMOR_FILE)) {

    if (
      REUSE_EXISTING_HEAVY_TABLES &&
      "tumor_meningioma_events" %in% dbListTables(con)
    ) {
      cat("\nReusing tumor_meningioma_events.\n")
    } else {

      cat("\nScanning tumor.csv for ICD-O meningioma...\n")
      t0 <- Sys.time()

      tumor_sql <- paste0(
  "
  CREATE OR REPLACE TABLE tumor_meningioma_events AS
  SELECT DISTINCT
    t.patient_id,
    TRY_CAST(STRPTIME(t.diagnosis_date, '%Y%m%d') AS DATE)
      AS diagnosis_date,
    UPPER(REPLACE(TRIM(t.tumor_site_code), '.', '')) AS site_norm,
    UPPER(TRIM(t.morphology_code)) AS morphology_code,
    t.source_id
  FROM read_csv_auto(
    '", sql_path(TUMOR_FILE), "',
    header=true,
    all_varchar=true,
    ignore_errors=false
  ) t
  JOIN analysis_cohort_primary_revised a
    ON a.patient_id=t.patient_id
  WHERE UPPER(REPLACE(TRIM(t.tumor_site_code), '.', ''))='C700'
    AND UPPER(TRIM(t.morphology_code)) LIKE '953%';
  "
      )

      dbExecute(con, tumor_sql)

      cat(
        "tumor.csv scan completed in ",
        round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1),
        " minutes.\n",
        sep=""
      )
    }

    dbExecute(con, "
    CREATE OR REPLACE TABLE analysis_cohort_registry_meningioma AS
    WITH o AS (
      SELECT
        a.patient_id,

        MIN(t.diagnosis_date) FILTER (
          WHERE t.diagnosis_date < a.index_date + INTERVAL 365 DAY
        ) AS registry_pre_landmark,

        MIN(t.diagnosis_date) FILTER (
          WHERE t.diagnosis_date >= a.index_date + INTERVAL 365 DAY
            AND t.diagnosis_date <= a.last_encounter_date
        ) AS registry_event_date,

        MIN(t.diagnosis_date) FILTER (
          WHERE UPPER(COALESCE(t.source_id,''))='CANCER REGISTRY'
            AND t.diagnosis_date < a.index_date + INTERVAL 365 DAY
        ) AS cancer_registry_pre_landmark,

        MIN(t.diagnosis_date) FILTER (
          WHERE UPPER(COALESCE(t.source_id,''))='CANCER REGISTRY'
            AND t.diagnosis_date >= a.index_date + INTERVAL 365 DAY
            AND t.diagnosis_date <= a.last_encounter_date
        ) AS cancer_registry_event_date

      FROM analysis_cohort_primary_revised a
      LEFT JOIN tumor_meningioma_events t
        ON t.patient_id=a.patient_id
      GROUP BY a.patient_id
    )
    SELECT
      a.*,

      CASE
        WHEN o.registry_event_date IS NOT NULL THEN 1 ELSE 0
      END AS registry_meningioma_event,

      o.registry_event_date,

      DATE_DIFF(
        'day',
        a.index_date + INTERVAL 365 DAY,
        CASE
          WHEN o.registry_event_date IS NOT NULL
          THEN o.registry_event_date
          ELSE a.last_encounter_date
        END
      ) AS registry_followup_days,

      CASE
        WHEN o.cancer_registry_event_date IS NOT NULL THEN 1 ELSE 0
      END AS cancer_registry_meningioma_event,

      o.cancer_registry_event_date,

      DATE_DIFF(
        'day',
        a.index_date + INTERVAL 365 DAY,
        CASE
          WHEN o.cancer_registry_event_date IS NOT NULL
          THEN o.cancer_registry_event_date
          ELSE a.last_encounter_date
        END
      ) AS cancer_registry_followup_days

    FROM analysis_cohort_primary_revised a
    JOIN o USING (patient_id)

    WHERE o.registry_pre_landmark IS NULL;
    ")
  }

  # ------------------------------------------------------------
  # 2. CLINICALLY TREATED MENINGIOMA VIA PROCEDURES
  # ------------------------------------------------------------
  # Meningioma-specific open resection:
  #   CPT 61512, 61519
  #
  # Cranial SRS treatment/management:
  #   CPT 61796-61799, 77371, 77372, 77432
  #   HCPCS G0173, G0251, G0339, G0340, G0563
  #
  # We require the treatment code from diagnosis through 180 days afterward.
  if (file.exists(PROCEDURE_FILE)) {

    if (
      REUSE_EXISTING_HEAVY_TABLES &&
      "meningioma_treatment_procedures" %in% dbListTables(con)
    ) {
      cat("\nReusing meningioma_treatment_procedures.\n")
    } else {

      cat("\nScanning procedure.csv for meningioma resection / cranial SRS...\n")
      t0 <- Sys.time()

      proc_sql <- paste0(
  "
  CREATE OR REPLACE TABLE meningioma_treatment_procedures AS
  SELECT DISTINCT
    p.patient_id,
    TRY_CAST(STRPTIME(p.date, '%Y%m%d') AS DATE) AS procedure_date,
    UPPER(TRIM(p.code_system)) AS code_system,
    UPPER(TRIM(p.code)) AS code,

    CASE
      WHEN UPPER(TRIM(p.code_system))='CPT'
       AND UPPER(TRIM(p.code)) IN ('61512','61519')
        THEN 'MENINGIOMA_RESECTION'

      WHEN (
        UPPER(TRIM(p.code_system))='CPT'
        AND UPPER(TRIM(p.code)) IN (
          '61796','61797','61798','61799',
          '77371','77372','77432'
        )
      )
      OR (
        UPPER(TRIM(p.code_system))='HCPCS'
        AND UPPER(TRIM(p.code)) IN (
          'G0173','G0251','G0339','G0340','G0563'
        )
      )
        THEN 'CRANIAL_SRS'

      ELSE NULL
    END AS treatment_type,

    p.source_id

  FROM read_csv_auto(
    '", sql_path(PROCEDURE_FILE), "',
    header=true,
    all_varchar=true,
    ignore_errors=false
  ) p
  JOIN analysis_cohort_primary_revised a
    ON a.patient_id=p.patient_id

  WHERE
    (
      UPPER(TRIM(p.code_system))='CPT'
      AND UPPER(TRIM(p.code)) IN (
        '61512','61519',
        '61796','61797','61798','61799',
        '77371','77372','77432'
      )
    )
    OR
    (
      UPPER(TRIM(p.code_system))='HCPCS'
      AND UPPER(TRIM(p.code)) IN (
        'G0173','G0251','G0339','G0340','G0563'
      )
    );
  "
      )

      dbExecute(con, proc_sql)

      cat(
        "procedure.csv scan completed in ",
        round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1),
        " minutes.\n",
        sep=""
      )
    }

    dbExecute(con, "
    CREATE OR REPLACE TABLE analysis_cohort_treated_meningioma AS
    WITH first_dx AS (
      SELECT
        a.patient_id,
        MIN(m.dx_date) FILTER (
          WHERE m.strict_meningioma=1
            AND m.dx_date >= a.index_date + INTERVAL 365 DAY
            AND m.dx_date <= a.last_encounter_date
        ) AS meningioma_dx_date
      FROM analysis_cohort_primary_revised a
      LEFT JOIN meningioma_dx_events_revised m
        ON m.patient_id=a.patient_id
      GROUP BY a.patient_id
    ),
    tx AS (
      SELECT
        d.patient_id,
        d.meningioma_dx_date,
        MIN(p.procedure_date) AS treatment_date,
        STRING_AGG(
          DISTINCT p.treatment_type,
          ', ' ORDER BY p.treatment_type
        ) AS treatment_types
      FROM first_dx d
      JOIN meningioma_treatment_procedures p
        ON p.patient_id=d.patient_id
       AND d.meningioma_dx_date IS NOT NULL
       AND p.procedure_date >= d.meningioma_dx_date
       AND p.procedure_date <= d.meningioma_dx_date + INTERVAL 180 DAY
      GROUP BY d.patient_id, d.meningioma_dx_date
    )
    SELECT
      a.*,
      d.meningioma_dx_date,
      tx.treatment_date,
      tx.treatment_types,

      CASE WHEN tx.treatment_date IS NOT NULL THEN 1 ELSE 0 END
        AS treated_meningioma_event,

      DATE_DIFF(
        'day',
        a.index_date + INTERVAL 365 DAY,
        CASE
          WHEN tx.treatment_date IS NOT NULL
          THEN d.meningioma_dx_date
          ELSE a.last_encounter_date
        END
      ) AS treated_meningioma_followup_days

    FROM analysis_cohort_primary_revised a
    LEFT JOIN first_dx d USING (patient_id)
    LEFT JOIN tx USING (patient_id);
    ")
  }

  # ------------------------------------------------------------
  # 3. ONCOLOGY-TREATMENT-TABLE ENDPOINT
  # ------------------------------------------------------------
  if (file.exists(ONCOLOGY_TREATMENT_FILE)) {

    if (
      REUSE_EXISTING_HEAVY_TABLES &&
      "oncology_treated_meningioma_events" %in% dbListTables(con)
    ) {
      cat("\nReusing oncology_treated_meningioma_events.\n")
    } else {

      cat("\nScanning oncology_treatment.csv for ICD-O meningioma treatment...\n")
      t0 <- Sys.time()

      onctx_sql <- paste0(
  "
  CREATE OR REPLACE TABLE oncology_treated_meningioma_events AS
  SELECT DISTINCT
    o.patient_id,
    TRY_CAST(STRPTIME(o.diagnosis_date, '%Y%m%d') AS DATE)
      AS diagnosis_date,
    TRY_CAST(STRPTIME(o.oncology_treatment_start_date, '%Y%m%d') AS DATE)
      AS treatment_start_date,
    UPPER(REPLACE(TRIM(o.tumor_site_code), '.', '')) AS site_norm,
    UPPER(TRIM(o.morphology_code)) AS morphology_code,
    o.oncology_treatment_code,
    t.code_description AS treatment_description,
    o.source_id
  FROM read_csv_auto(
    '", sql_path(ONCOLOGY_TREATMENT_FILE), "',
    header=true,
    all_varchar=true,
    ignore_errors=false
  ) o
  LEFT JOIN standardized_terminology t
    ON t.code_system=o.oncology_treatment_code_system
   AND t.code=o.oncology_treatment_code
  JOIN analysis_cohort_primary_revised a
    ON a.patient_id=o.patient_id
  WHERE UPPER(REPLACE(TRIM(o.tumor_site_code), '.', ''))='C700'
    AND UPPER(TRIM(o.morphology_code)) LIKE '953%'
    AND o.oncology_treatment_start_date IS NOT NULL
    AND o.oncology_treatment_start_date <> '';
  "
      )

      dbExecute(con, onctx_sql)

      cat(
        "oncology_treatment.csv scan completed in ",
        round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1),
        " minutes.\n",
        sep=""
      )
    }

    dbExecute(con, "
    CREATE OR REPLACE TABLE analysis_cohort_oncology_treated_registry AS
    WITH o AS (
      SELECT
        a.patient_id,

        MIN(x.diagnosis_date) FILTER (
          WHERE x.diagnosis_date < a.index_date + INTERVAL 365 DAY
        ) AS pre_landmark_dx,

        MIN(x.diagnosis_date) FILTER (
          WHERE x.diagnosis_date >= a.index_date + INTERVAL 365 DAY
            AND x.diagnosis_date <= a.last_encounter_date
        ) AS treated_registry_event_date

      FROM analysis_cohort_primary_revised a
      LEFT JOIN oncology_treated_meningioma_events x
        ON x.patient_id=a.patient_id
      GROUP BY a.patient_id
    )
    SELECT
      a.*,

      CASE
        WHEN o.treated_registry_event_date IS NOT NULL THEN 1 ELSE 0
      END AS oncology_treated_registry_event,

      o.treated_registry_event_date,

      DATE_DIFF(
        'day',
        a.index_date + INTERVAL 365 DAY,
        CASE
          WHEN o.treated_registry_event_date IS NOT NULL
          THEN o.treated_registry_event_date
          ELSE a.last_encounter_date
        END
      ) AS oncology_treated_registry_followup_days

    FROM analysis_cohort_primary_revised a
    JOIN o USING (patient_id)
    WHERE o.pre_landmark_dx IS NULL;
    ")
  }

  # ------------------------------------------------------------
  # 4. MEASURED BMI
  # ------------------------------------------------------------
  if (RUN_BMI_SCAN && file.exists(VITALS_FILE)) {

    bmi_codes <- as.data.table(dbGetQuery(con, "
    SELECT
      code_system,
      code,
      code_description,
      unit
    FROM standardized_terminology
    WHERE UPPER(code_system)='LOINC'
      AND (
        code='39156-5'
        OR LOWER(code_description) LIKE '%body mass index%'
      )
    ORDER BY code;
    "))

    fwrite(bmi_codes, file.path(QC_DIR, "09_bmi_loinc_codes.csv"))

    if (nrow(bmi_codes) == 0L) {
      warning("No BMI LOINC terms found in standardized_terminology; BMI scan skipped.")
    } else {

      code_sql <- paste(
        sprintf("'%s'", gsub("'", "''", unique(bmi_codes$code), fixed=TRUE)),
        collapse=","
      )

      if (
        REUSE_EXISTING_HEAVY_TABLES &&
        "primary_preindex_bmi" %in% dbListTables(con)
      ) {
        cat("\nReusing primary_preindex_bmi.\n")
      } else {

        cat("\n============================================================\n")
        cat("Scanning vitals_signs.csv for measured BMI...\n")
        cat("============================================================\n")
        t0 <- Sys.time()

        bmi_sql <- paste0(
  "
  CREATE OR REPLACE TABLE primary_preindex_bmi AS
  WITH v AS (
    SELECT
      patient_id,
      TRY_CAST(STRPTIME(date, '%Y%m%d') AS DATE) AS vital_date,
      TRY_CAST(value AS DOUBLE) AS bmi_value,
      code,
      units_of_measure,
      source_id
    FROM read_csv_auto(
      '", sql_path(VITALS_FILE), "',
      header=true,
      all_varchar=true,
      ignore_errors=false
    )
    WHERE UPPER(TRIM(code_system))='LOINC'
      AND code IN (", code_sql, ")
  ),
  eligible AS (
    SELECT
      a.patient_id,
      v.vital_date,
      v.bmi_value,
      v.code,
      v.units_of_measure,
      v.source_id,
      ROW_NUMBER() OVER (
        PARTITION BY a.patient_id
        ORDER BY v.vital_date DESC
      ) AS rn
    FROM analysis_cohort_primary_revised a
    JOIN v
      ON v.patient_id=a.patient_id
     AND v.vital_date >= a.index_date - INTERVAL 365 DAY
     AND v.vital_date < a.index_date
    WHERE v.bmi_value BETWEEN 10 AND 80
  )
  SELECT
    patient_id,
    vital_date AS bmi_date,
    bmi_value AS bmi_preindex,
    code AS bmi_code,
    units_of_measure AS bmi_units,
    source_id AS bmi_source_id
  FROM eligible
  WHERE rn=1;
  "
        )

        dbExecute(con, bmi_sql)

        cat(
          "BMI scan completed in ",
          round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1),
          " minutes.\n",
          sep=""
        )
      }

      dbExecute(con, "
      CREATE OR REPLACE TABLE analysis_cohort_primary_bmi AS
      SELECT
        a.*,
        b.bmi_date,
        b.bmi_preindex,
        CASE WHEN b.bmi_preindex IS NULL THEN 1 ELSE 0 END AS bmi_missing
      FROM analysis_cohort_primary_revised a
      LEFT JOIN primary_preindex_bmi b USING (patient_id);
      ")

      bmi_qc <- as.data.table(dbGetQuery(con, "
      SELECT
        exposure,
        COUNT(*) AS patients,
        SUM(CASE WHEN bmi_preindex IS NOT NULL THEN 1 ELSE 0 END)
          AS bmi_observed,
        ROUND(
          100.0 * SUM(CASE WHEN bmi_preindex IS NOT NULL THEN 1 ELSE 0 END)
          / COUNT(*),
          1
        ) AS bmi_observed_pct,
        MEDIAN(bmi_preindex) AS median_bmi
      FROM analysis_cohort_primary_bmi
      GROUP BY exposure
      ORDER BY exposure;
      "))

      fwrite(bmi_qc, file.path(QC_DIR, "09_bmi_completeness.csv"))
    }
  }

  # ------------------------------------------------------------
  # 5. QC
  # ------------------------------------------------------------
  qc_list <- list()

  if ("analysis_cohort_registry_meningioma" %in% dbListTables(con)) {
    qc_list[["registry"]] <- as.data.table(dbGetQuery(con, "
    SELECT
      exposure,
      COUNT(*) AS patients,
      SUM(registry_meningioma_event) AS tumor_table_events,
      SUM(cancer_registry_meningioma_event) AS cancer_registry_events
    FROM analysis_cohort_registry_meningioma
    GROUP BY exposure
    ORDER BY exposure;
    "))
    fwrite(
      qc_list[["registry"]],
      file.path(QC_DIR, "09_registry_outcome_counts.csv")
    )
  }

  if ("analysis_cohort_treated_meningioma" %in% dbListTables(con)) {
    qc_list[["treated"]] <- as.data.table(dbGetQuery(con, "
    SELECT
      exposure,
      COUNT(*) AS patients,
      SUM(treated_meningioma_event) AS treated_events
    FROM analysis_cohort_treated_meningioma
    GROUP BY exposure
    ORDER BY exposure;
    "))
    fwrite(
      qc_list[["treated"]],
      file.path(QC_DIR, "09_treated_outcome_counts.csv")
    )

    tx_types <- as.data.table(dbGetQuery(con, "
    SELECT
      treatment_types,
      COUNT(*) AS patients
    FROM analysis_cohort_treated_meningioma
    WHERE treated_meningioma_event=1
    GROUP BY treatment_types
    ORDER BY patients DESC;
    "))
    fwrite(
      tx_types,
      file.path(QC_DIR, "09_treated_outcome_types.csv")
    )
  }

  if ("analysis_cohort_oncology_treated_registry" %in% dbListTables(con)) {
    qc_list[["oncology_treated"]] <- as.data.table(dbGetQuery(con, "
    SELECT
      exposure,
      COUNT(*) AS patients,
      SUM(oncology_treated_registry_event) AS events
    FROM analysis_cohort_oncology_treated_registry
    GROUP BY exposure
    ORDER BY exposure;
    "))
    fwrite(
      qc_list[["oncology_treated"]],
      file.path(QC_DIR, "09_oncology_treated_registry_counts.csv")
    )
  }

  cat("\n================ SECONDARY ENDPOINT QC ==============\n")
  for (nm in names(qc_list)) {
    cat("\n--- ", nm, " ---\n", sep="")
    print(qc_list[[nm]])
  }

  cat("\nStage 09 complete.\n")

  dbDisconnect(con, shutdown = TRUE)
}

run_stage_09()
