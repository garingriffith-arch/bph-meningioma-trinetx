# NOTE: This stage is wrapped in a function so on.exit() is tied to
# the stage function rather than source()/eval(). This prevents the
# DuckDB connection from being closed immediately when sourced.

run_stage_07 <- function() {
  # ============================================================
  # 07_build_core_final_sensitivities.R
  #
  # Core literature-informed sensitivity cohorts that reuse the
  # existing DuckDB work wherever possible.
  #
  # Builds:
  #   1) 1-year sustained-treatment landmark cohort
  #   2) stricter 1-year sustained-treatment cohort
  #   3) 3-year persistent-treatment landmark cohort
  #   4) high-specificity intracranial meningioma outcomes
  #   5) latency-analysis ingredients
  #   6) death competing-risk ingredients
  #   7) medication_drug days-supply audit / exploratory cohort
  #
  # Heavy new scan:
  #   medication_drug.csv ONCE (checkpointed)
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

  MED_DRUG_FILE <- file.path(RAW_DIR, "medication_drug.csv")

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

  needed <- c(
    "analysis_cohort_primary_revised",
    "revised_eligible_base",
    "study_med_events",
    "study_med_events_raw",
    "meningioma_dx_events_revised",
    "patient_core"
  )
  missing <- setdiff(needed, dbListTables(con))
  if (length(missing) > 0L) {
    stop(
      "Missing required prior-stage tables: ",
      paste(missing, collapse = ", ")
    )
  }

  # ------------------------------------------------------------
  # 1. MEDICATION PERSISTENCE METRICS
  # ------------------------------------------------------------
  dbExecute(con, "
  CREATE OR REPLACE TABLE final_med_persistence AS
  WITH b AS (
    SELECT
      patient_id,
      exposure,
      index_date,
      last_encounter_date
    FROM revised_eligible_base
  ),
  j AS (
    SELECT
      b.patient_id,
      b.exposure,
      b.index_date,
      b.last_encounter_date,
      e.drug,
      e.drug_class,
      e.start_date,

      CASE
        WHEN b.exposure='5ARI' AND e.drug_class='5ARI' THEN 1
        WHEN b.exposure='URO_ALPHA' AND e.drug_class='URO_ALPHA' THEN 1
        ELSE 0
      END AS same_class,

      CASE
        WHEN b.exposure='5ARI' AND e.drug_class='URO_ALPHA' THEN 1
        WHEN b.exposure='URO_ALPHA' AND e.drug_class='5ARI' THEN 1
        ELSE 0
      END AS opposite_class

    FROM b
    LEFT JOIN study_med_events e
      ON e.patient_id=b.patient_id
  )
  SELECT
    patient_id,

    COUNT(DISTINCT start_date) FILTER (
      WHERE same_class=1
        AND start_date > index_date
        AND start_date <= index_date + INTERVAL 365 DAY
    ) AS same_dates_0_365,

    COUNT(DISTINCT start_date) FILTER (
      WHERE same_class=1
        AND start_date >= index_date + INTERVAL 90 DAY
        AND start_date <= index_date + INTERVAL 365 DAY
    ) AS same_dates_90_365,

    COUNT(DISTINCT start_date) FILTER (
      WHERE same_class=1
        AND start_date >= index_date + INTERVAL 90 DAY
        AND start_date < index_date + INTERVAL 181 DAY
    ) AS same_dates_90_180,

    COUNT(DISTINCT start_date) FILTER (
      WHERE same_class=1
        AND start_date >= index_date + INTERVAL 181 DAY
        AND start_date <= index_date + INTERVAL 365 DAY
    ) AS same_dates_181_365,

    COUNT(DISTINCT start_date) FILTER (
      WHERE opposite_class=1
        AND start_date > index_date
        AND start_date <= index_date + INTERVAL 365 DAY
    ) AS opposite_dates_0_365,

    COUNT(DISTINCT start_date) FILTER (
      WHERE same_class=1
        AND start_date > index_date
        AND start_date <= index_date + INTERVAL 1095 DAY
    ) AS same_dates_0_1095,

    COUNT(DISTINCT start_date) FILTER (
      WHERE same_class=1
        AND start_date > index_date
        AND start_date <= index_date + INTERVAL 365 DAY
    ) AS same_dates_y1,

    COUNT(DISTINCT start_date) FILTER (
      WHERE same_class=1
        AND start_date > index_date + INTERVAL 365 DAY
        AND start_date <= index_date + INTERVAL 730 DAY
    ) AS same_dates_y2,

    COUNT(DISTINCT start_date) FILTER (
      WHERE same_class=1
        AND start_date > index_date + INTERVAL 730 DAY
        AND start_date <= index_date + INTERVAL 1095 DAY
    ) AS same_dates_y3,

    COUNT(DISTINCT start_date) FILTER (
      WHERE opposite_class=1
        AND start_date > index_date
        AND start_date <= index_date + INTERVAL 1095 DAY
    ) AS opposite_dates_0_1095

  FROM j
  GROUP BY patient_id;
  ")

  # ------------------------------------------------------------
  # 2. 1-YEAR SUSTAINED-TREATMENT LANDMARK
  # ------------------------------------------------------------
  # Primary sustained definition:
  #   - >=1 repeat same-class record between days 90-365
  #   - no opposite study class during days 1-365
  #   - already tumor-free during the 1-y lag by virtue of the revised primary
  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_sustained_1y AS
  SELECT
    a.*,
    p.same_dates_0_365,
    p.same_dates_90_365,
    p.same_dates_90_180,
    p.same_dates_181_365,
    p.opposite_dates_0_365
  FROM analysis_cohort_primary_revised a
  JOIN final_med_persistence p USING (patient_id)
  WHERE p.same_dates_90_365 >= 1
    AND p.opposite_dates_0_365 = 0;
  ")

  # Stricter version:
  # evidence in both halves of the post-index year.
  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_sustained_1y_strict AS
  SELECT *
  FROM analysis_cohort_sustained_1y
  WHERE same_dates_90_180 >= 1
    AND same_dates_181_365 >= 1;
  ")

  # ------------------------------------------------------------
  # 3. 3-YEAR PERSISTENT-TREATMENT LANDMARK
  # ------------------------------------------------------------
  # Require:
  #   - >=1 same-class medication record in EACH of years 1, 2, 3
  #   - no opposite class during first 3 years
  #   - no broad meningeal tumor before the 3-y landmark
  #   - positive observation after the landmark
  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_persistent_3y AS
  WITH tumor3 AS (
    SELECT
      b.patient_id,
      MIN(m.dx_date) FILTER (
        WHERE m.broad_meningeal_tumor=1
          AND m.dx_date >= b.index_date
          AND m.dx_date < b.index_date + INTERVAL 1095 DAY
      ) AS broad_tumor_pre_3y,

      MIN(m.dx_date) FILTER (
        WHERE m.strict_meningioma=1
          AND m.dx_date >= b.index_date + INTERVAL 1095 DAY
          AND m.dx_date <= b.last_encounter_date
      ) AS strict_event_post_3y

    FROM revised_eligible_base b
    LEFT JOIN meningioma_dx_events_revised m
      ON m.patient_id=b.patient_id
    GROUP BY b.patient_id
  )
  SELECT
    b.*,
    p.same_dates_y1,
    p.same_dates_y2,
    p.same_dates_y3,
    p.opposite_dates_0_1095,

    b.index_date + INTERVAL 1095 DAY AS landmark_date,

    CASE
      WHEN t.strict_event_post_3y IS NOT NULL THEN 1 ELSE 0
    END AS meningioma_event_strict,

    t.strict_event_post_3y AS incident_meningioma_date_strict,

    DATE_DIFF(
      'day',
      b.index_date + INTERVAL 1095 DAY,
      CASE
        WHEN t.strict_event_post_3y IS NOT NULL
        THEN t.strict_event_post_3y
        ELSE b.last_encounter_date
      END
    ) AS followup_days_strict

  FROM revised_eligible_base b
  JOIN final_med_persistence p USING (patient_id)
  JOIN tumor3 t USING (patient_id)

  WHERE p.same_dates_y1 >= 1
    AND p.same_dates_y2 >= 1
    AND p.same_dates_y3 >= 1
    AND p.opposite_dates_0_1095 = 0
    AND t.broad_tumor_pre_3y IS NULL
    AND b.last_encounter_date > b.index_date + INTERVAL 1095 DAY;
  ")

  # ------------------------------------------------------------
  # 4. HIGH-SPECIFICITY MENINGIOMA OUTCOMES
  # ------------------------------------------------------------
  # Specific intracranial meningioma:
  #   ICD-10-CM D32.0 + historical ICD-9-CM 225.2
  #
  # Pure ICD-10 sensitivity:
  #   D32.0 only, with index restricted to ICD-10 era.
  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_specific_meningioma AS
  WITH o AS (
    SELECT
      a.patient_id,

      MIN(m.dx_date) FILTER (
        WHERE (
              (m.code_system='ICD-10-CM' AND m.code_norm='D320')
           OR (m.code_system='ICD-9-CM'  AND m.code_norm='2252')
        )
        AND m.dx_date >= a.index_date + INTERVAL 365 DAY
        AND m.dx_date <= a.last_encounter_date
      ) AS specific_event_date,

      MIN(m.dx_date) FILTER (
        WHERE m.code_system='ICD-10-CM'
          AND m.code_norm='D320'
          AND m.dx_date >= a.index_date + INTERVAL 365 DAY
          AND m.dx_date <= a.last_encounter_date
      ) AS d320_event_date

    FROM analysis_cohort_primary_revised a
    LEFT JOIN meningioma_dx_events_revised m
      ON m.patient_id=a.patient_id
    GROUP BY a.patient_id
  )
  SELECT
    a.*,

    CASE WHEN o.specific_event_date IS NOT NULL THEN 1 ELSE 0 END
      AS meningioma_event_specific,

    o.specific_event_date,

    DATE_DIFF(
      'day',
      a.index_date + INTERVAL 365 DAY,
      CASE
        WHEN o.specific_event_date IS NOT NULL
        THEN o.specific_event_date
        ELSE a.last_encounter_date
      END
    ) AS followup_days_specific,

    CASE WHEN o.d320_event_date IS NOT NULL THEN 1 ELSE 0 END
      AS meningioma_event_d320,

    o.d320_event_date,

    DATE_DIFF(
      'day',
      a.index_date + INTERVAL 365 DAY,
      CASE
        WHEN o.d320_event_date IS NOT NULL
        THEN o.d320_event_date
        ELSE a.last_encounter_date
      END
    ) AS followup_days_d320

  FROM analysis_cohort_primary_revised a
  JOIN o USING (patient_id);
  ")

  # ------------------------------------------------------------
  # 5. LATENCY INGREDIENTS FOR CLEAN COHORT
  # ------------------------------------------------------------
  # This table is only an ingredient table. The less selected latency cohort
  # (including patients with <1 y follow-up) is built in Stage 08 after a fresh
  # diagnosis scan. This table provides a checkpoint / comparison using the
  # already-built revised eligible base.
  dbExecute(con, "
  CREATE OR REPLACE TABLE latency_revised_base_checkpoint AS
  WITH o AS (
    SELECT
      b.patient_id,
      MIN(m.dx_date) FILTER (
        WHERE m.strict_meningioma=1
          AND m.dx_date >= b.index_date
          AND m.dx_date <= b.last_encounter_date
      ) AS strict_event_post_index
    FROM revised_eligible_base b
    LEFT JOIN meningioma_dx_events_revised m
      ON m.patient_id=b.patient_id
    GROUP BY b.patient_id
  )
  SELECT
    b.*,
    CASE WHEN o.strict_event_post_index IS NOT NULL THEN 1 ELSE 0 END
      AS meningioma_event_strict,
    o.strict_event_post_index AS incident_meningioma_date_strict,
    DATE_DIFF(
      'day',
      b.index_date,
      CASE
        WHEN o.strict_event_post_index IS NOT NULL
        THEN o.strict_event_post_index
        ELSE b.last_encounter_date
      END
    ) AS followup_days_strict
  FROM revised_eligible_base b
  JOIN o USING (patient_id);
  ")

  # ------------------------------------------------------------
  # 6. COMPETING-RISK INGREDIENTS
  # ------------------------------------------------------------
  # month_year_death is intentionally coarse/deidentified. We use the first
  # day of the provided death month only as a marker for competing-risk
  # sensitivity analysis, and write validity QC.
  dbExecute(con, "
  CREATE OR REPLACE TABLE patient_death_marker AS
  SELECT
    patient_id,
    month_year_death,

    CASE
      WHEN month_year_death IS NOT NULL
       AND CAST(FLOOR(month_year_death / 100) AS INTEGER)
             BETWEEN 1900 AND 2100
       AND CAST(month_year_death % 100 AS INTEGER)
             BETWEEN 1 AND 12
      THEN MAKE_DATE(
        CAST(FLOOR(month_year_death / 100) AS INTEGER),
        CAST(month_year_death % 100 AS INTEGER),
        1
      )
      ELSE NULL
    END AS death_month_marker

  FROM patient_core;
  ")

  death_qc <- as.data.table(dbGetQuery(con, "
  SELECT
    COUNT(*) AS patient_rows,
    SUM(CASE WHEN month_year_death IS NOT NULL THEN 1 ELSE 0 END)
      AS nonmissing_death_month,
    SUM(CASE WHEN month_year_death IS NOT NULL
              AND death_month_marker IS NULL THEN 1 ELSE 0 END)
      AS invalid_death_month,
    SUM(CASE WHEN death_month_marker IS NOT NULL THEN 1 ELSE 0 END)
      AS valid_death_month
  FROM patient_death_marker;
  "))
  fwrite(death_qc, file.path(QC_DIR, "07_death_marker_qc.csv"))

  dbExecute(con, "
  CREATE OR REPLACE TABLE analysis_cohort_competing_risk AS
  SELECT
    a.*,
    d.death_month_marker,

    CASE
      WHEN a.incident_meningioma_date_strict IS NOT NULL
       AND (
         d.death_month_marker IS NULL
         OR a.incident_meningioma_date_strict <= d.death_month_marker
       )
       AND a.incident_meningioma_date_strict <= a.last_encounter_date
        THEN 1

      WHEN d.death_month_marker IS NOT NULL
       AND d.death_month_marker > a.index_date + INTERVAL 365 DAY
       AND d.death_month_marker < a.last_encounter_date
        THEN 2

      ELSE 0
    END AS competing_status,

    DATE_DIFF(
      'day',
      a.index_date + INTERVAL 365 DAY,
      CASE
        WHEN a.incident_meningioma_date_strict IS NOT NULL
         AND (
           d.death_month_marker IS NULL
           OR a.incident_meningioma_date_strict <= d.death_month_marker
         )
         AND a.incident_meningioma_date_strict <= a.last_encounter_date
          THEN a.incident_meningioma_date_strict

        WHEN d.death_month_marker IS NOT NULL
         AND d.death_month_marker > a.index_date + INTERVAL 365 DAY
         AND d.death_month_marker < a.last_encounter_date
          THEN d.death_month_marker

        ELSE a.last_encounter_date
      END
    ) AS competing_time_days

  FROM analysis_cohort_primary_revised a
  LEFT JOIN patient_death_marker d USING (patient_id);
  ")

  # ------------------------------------------------------------
  # 7. MEDICATION_DRUG / DAYS-SUPPLY AUDIT
  # ------------------------------------------------------------
  if (!file.exists(MED_DRUG_FILE)) {
    warning("medication_drug.csv not found. Days-supply audit skipped.")
  } else if (
    REUSE_EXISTING_HEAVY_TABLES &&
    "study_med_days_supply" %in% dbListTables(con)
  ) {
    cat("\nReusing existing study_med_days_supply table.\n")
  } else {

    cat("\n============================================================\n")
    cat("Scanning medication_drug.csv once for days_supply...\n")
    cat("============================================================\n")
    t0 <- Sys.time()

    # study_med_events_raw retained unique_id from Stage 1, so we can recover
    # product-level days_supply without rescanning medication_ingredient.csv.
    days_sql <- paste0(
  "
  CREATE OR REPLACE TABLE study_med_days_supply AS
  WITH study_ids AS (
    SELECT DISTINCT
      patient_id,
      unique_id,
      drug,
      drug_class,
      start_date
    FROM study_med_events_raw
    WHERE unique_id IS NOT NULL
      AND unique_id <> ''
  ),
  md AS (
    SELECT
      patient_id,
      unique_id,
      TRY_CAST(STRPTIME(start_date, '%Y%m%d') AS DATE) AS drug_start_date,
      days_supply,
      quantity_dispensed,
      source_id,

      TRY_CAST(
        NULLIF(
          REGEXP_EXTRACT(TRIM(days_supply), '([0-9]+[.]?[0-9]*)', 1),
          ''
        ) AS DOUBLE
      ) AS days_supply_num

    FROM read_csv_auto(
      '", sql_path(MED_DRUG_FILE), "',
      header=true,
      all_varchar=true,
      ignore_errors=false
    )
  )
  SELECT
    s.patient_id,
    s.unique_id,
    s.drug,
    s.drug_class,
    s.start_date,

    MAX(md.days_supply_num) AS days_supply_num,
    MIN(md.days_supply) AS days_supply_raw,
    MIN(md.quantity_dispensed) AS quantity_dispensed,
    MIN(md.source_id) AS source_id,

    COUNT(*) AS matched_drug_rows

  FROM study_ids s
  JOIN md
    ON md.patient_id=s.patient_id
   AND md.unique_id=s.unique_id

  GROUP BY
    s.patient_id,
    s.unique_id,
    s.drug,
    s.drug_class,
    s.start_date;
  "
    )

    dbExecute(con, days_sql)

    try(dbExecute(
      con,
      "CREATE INDEX IF NOT EXISTS idx_study_days_patient_date
       ON study_med_days_supply(patient_id, start_date);"
    ), silent = TRUE)

    cat(
      "medication_drug.csv scan completed in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      " minutes.\n",
      sep = ""
    )
  }

  if ("study_med_days_supply" %in% dbListTables(con)) {

    ds_qc <- as.data.table(dbGetQuery(con, "
    SELECT
      drug,
      source_id,
      COUNT(*) AS matched_events,
      SUM(CASE WHEN days_supply_num IS NOT NULL
                AND days_supply_num > 0
                AND days_supply_num <= 365
               THEN 1 ELSE 0 END) AS valid_days_supply_events,
      ROUND(
        100.0 *
        SUM(CASE WHEN days_supply_num IS NOT NULL
                  AND days_supply_num > 0
                  AND days_supply_num <= 365
                 THEN 1 ELSE 0 END)
        / COUNT(*),
        1
      ) AS valid_days_supply_pct,
      MEDIAN(
        CASE WHEN days_supply_num > 0 AND days_supply_num <= 365
             THEN days_supply_num END
      ) AS median_valid_days_supply
    FROM study_med_days_supply
    GROUP BY drug, source_id
    ORDER BY drug, matched_events DESC;
    "))

    fwrite(ds_qc, file.path(QC_DIR, "07_days_supply_qc.csv"))

    # Patient-level first-year days-supply proxy.
    # This is explicitly exploratory: overlapping prescriptions are not
    # de-overlapped, so the sum is capped at 365 rather than called true PDC.
    dbExecute(con, "
    CREATE OR REPLACE TABLE primary_days_supply_1y AS
    WITH j AS (
      SELECT
        a.patient_id,
        a.exposure,
        a.index_date,
        d.drug_class,
        d.start_date,
        d.days_supply_num,

        CASE
          WHEN a.exposure='5ARI' AND d.drug_class='5ARI' THEN 1
          WHEN a.exposure='URO_ALPHA' AND d.drug_class='URO_ALPHA' THEN 1
          ELSE 0
        END AS same_class,

        CASE
          WHEN a.exposure='5ARI' AND d.drug_class='URO_ALPHA' THEN 1
          WHEN a.exposure='URO_ALPHA' AND d.drug_class='5ARI' THEN 1
          ELSE 0
        END AS opposite_class

      FROM analysis_cohort_primary_revised a
      LEFT JOIN study_med_days_supply d
        ON d.patient_id=a.patient_id
       AND d.start_date >= a.index_date
       AND d.start_date <= a.index_date + INTERVAL 365 DAY
    )
    SELECT
      patient_id,

      COUNT(*) FILTER (WHERE same_class=1) AS same_events_365,

      COUNT(*) FILTER (
        WHERE same_class=1
          AND days_supply_num IS NOT NULL
          AND days_supply_num > 0
          AND days_supply_num <= 365
      ) AS same_valid_ds_events_365,

      LEAST(
        365.0,
        COALESCE(
          SUM(
            CASE
              WHEN same_class=1
               AND days_supply_num > 0
               AND days_supply_num <= 365
              THEN days_supply_num
              ELSE 0
            END
          ),
          0
        )
      ) AS approx_days_supply_sum_365,

      COUNT(*) FILTER (WHERE opposite_class=1) AS opposite_events_365

    FROM j
    GROUP BY patient_id;
    ")

    dbExecute(con, "
    CREATE OR REPLACE TABLE analysis_cohort_days_supply_exploratory AS
    SELECT
      a.*,
      d.same_events_365,
      d.same_valid_ds_events_365,
      d.approx_days_supply_sum_365,
      d.opposite_events_365,

      CASE
        WHEN d.same_events_365 > 0
         AND d.same_valid_ds_events_365 = d.same_events_365
        THEN 1 ELSE 0
      END AS complete_days_supply_first_year

    FROM analysis_cohort_primary_revised a
    LEFT JOIN primary_days_supply_1y d USING (patient_id)

    WHERE d.same_events_365 > 0
      AND d.same_valid_ds_events_365 = d.same_events_365
      AND d.approx_days_supply_sum_365 >= 180
      AND d.opposite_events_365 = 0;
    ")

    ds_patient_qc <- as.data.table(dbGetQuery(con, "
    SELECT
      a.exposure,
      COUNT(*) AS primary_patients,
      SUM(CASE WHEN d.same_events_365 > 0 THEN 1 ELSE 0 END)
        AS patients_with_matched_same_class_drug_record,
      SUM(CASE WHEN d.same_events_365 > 0
                AND d.same_valid_ds_events_365 = d.same_events_365
               THEN 1 ELSE 0 END)
        AS patients_complete_days_supply_first_year,
      SUM(CASE WHEN d.same_events_365 > 0
                AND d.same_valid_ds_events_365 = d.same_events_365
                AND d.approx_days_supply_sum_365 >= 180
                AND d.opposite_events_365 = 0
               THEN 1 ELSE 0 END)
        AS patients_exploratory_sustained_ds
    FROM analysis_cohort_primary_revised a
    LEFT JOIN primary_days_supply_1y d USING (patient_id)
    GROUP BY a.exposure
    ORDER BY a.exposure;
    "))

    fwrite(
      ds_patient_qc,
      file.path(QC_DIR, "07_days_supply_patient_qc.csv")
    )
  }

  # ------------------------------------------------------------
  # 8. QC COUNTS
  # ------------------------------------------------------------
  count_table <- function(tbl, event_col = "meningioma_event_strict") {
    sql <- paste0(
      "SELECT exposure, COUNT(*) AS patients, ",
      if (event_col %in% dbListFields(con, tbl))
        paste0("SUM(", event_col, ") AS events")
      else
        "NULL AS events",
      " FROM ", tbl,
      " GROUP BY exposure ORDER BY exposure;"
    )
    as.data.table(dbGetQuery(con, sql))
  }

  qc_sustained <- count_table("analysis_cohort_sustained_1y")
  qc_sustained_strict <- count_table("analysis_cohort_sustained_1y_strict")
  qc_persistent3 <- count_table("analysis_cohort_persistent_3y")
  qc_specific <- as.data.table(dbGetQuery(con, "
  SELECT
    exposure,
    COUNT(*) AS patients,
    SUM(meningioma_event_specific) AS specific_events,
    SUM(meningioma_event_d320) AS d320_events
  FROM analysis_cohort_specific_meningioma
  GROUP BY exposure
  ORDER BY exposure;
  "))
  qc_competing <- as.data.table(dbGetQuery(con, "
  SELECT
    exposure,
    competing_status,
    COUNT(*) AS patients
  FROM analysis_cohort_competing_risk
  GROUP BY exposure, competing_status
  ORDER BY exposure, competing_status;
  "))

  fwrite(qc_sustained, file.path(QC_DIR, "07_sustained_1y_counts.csv"))
  fwrite(
    qc_sustained_strict,
    file.path(QC_DIR, "07_sustained_1y_strict_counts.csv")
  )
  fwrite(qc_persistent3, file.path(QC_DIR, "07_persistent_3y_counts.csv"))
  fwrite(qc_specific, file.path(QC_DIR, "07_specific_outcome_counts.csv"))
  fwrite(qc_competing, file.path(QC_DIR, "07_competing_risk_counts.csv"))

  cat("\n================ 1-Y SUSTAINED =====================\n")
  print(qc_sustained)

  cat("\n================ STRICT 1-Y SUSTAINED ==============\n")
  print(qc_sustained_strict)

  cat("\n================ 3-Y PERSISTENT ====================\n")
  print(qc_persistent3)

  cat("\n================ SPECIFIC OUTCOMES ==================\n")
  print(qc_specific)

  cat("\n================ COMPETING-RISK STATUS =============\n")
  print(qc_competing)

  cat("\nStage 07 complete.\n")

  dbDisconnect(con, shutdown = TRUE)
}

run_stage_07()
