# ============================================================
# 17_make_final_manuscript_outputs.R
#
# Finasteride / Dutasteride and Incident Meningioma
# FINAL manuscript-ready tables + figures
# PATCH v9 (pre-resident): improved borderline FDR display, source-neutral notes, balance/latency clarifications.
# Inherits v8 visual QA and boundary-aligned 5-y analysis.
#
# Incorporates finalized analyses through Stage 16:
#   - canonical boundary-fixed primary cohort
#   - target-trial-aligned analyses
#   - EHR-continuity IPCW
#   - internal phenotype corroboration
#   - prespecified falsification outcomes
#   - expanded exploratory BH-FDR family
#   - high-specificity persistence analyses
#
# Main outputs:
#   Table 1. Baseline characteristics
#   Table 2. Primary and robustness analyses
#   Table 3. Exploratory sustained/persistent exposure analyses
#
# Supplement:
#   Table S1. Target-trial specification
#   Table S2. Covariate balance
#   Table S3. Latency analyses
#   Table S4. Internal phenotype corroboration
#   Table S5. Prespecified falsification outcomes
#   Table S6. Complete exploratory multiple-comparison family
#   Table S7. Rare / clinically consequential outcome definitions
#
# Figures:
#   Figure 1. Cohort flow
#   Figure 2. Primary/robustness forest plot
#   Figure 3. Exploratory persistence forest plot
#   Figure S1. Love plot
#
# Exports:
#   - editable CSVs
#   - combined Word document
#   - 600-dpi PNG/TIFF + vector PDF figures
#
# Typography:
#   Word table body: Arial 10
#   Word table titles: Times New Roman 12
#   Graphics: generic sans (Windows/PDF safe)
#
# IMPORTANT:
# The finalized inferential results are read from Stages 06/10/12/15/16.
# This script deliberately refits ONLY three already-defined sensitivity
# models on the canonical boundary-fixed primary cohort:
#   (1) high-specificity D32.0/ICD-9 225.2,
#   (2) index 2012 onward,
#   (3) ICD-10-era index.
# This is solely to remove the one exact-landmark D32.9 patient from those
# denominators and keep manuscript counts internally consistent.
# No new hypothesis or subgroup is introduced here.
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------
# 0. PATHS
# ------------------------------------------------------------
R_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"
WORK_DIR <- "C:/TriNetX/bph_m_work"
DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")

ANALYSIS_DIR <- file.path(WORK_DIR, "analysis")
FINAL_DIR <- file.path(WORK_DIR, "analysis_final")
SIGNAL_DIR <- file.path(WORK_DIR, "analysis_signal_interrogation")
METHOD_DIR <- file.path(WORK_DIR, "analysis_methodologic_validity")
METHOD_QC_DIR <- file.path(WORK_DIR, "qc_methodologic_validity")
STAGE16_DIR <- file.path(WORK_DIR, "analysis_stage16_high_specificity")
STAGE16_QC_DIR <- file.path(WORK_DIR, "qc_stage16_high_specificity")
QC_DIR <- file.path(WORK_DIR, "qc")

OUT_DIR <- file.path(WORK_DIR, "manuscript_ready_final")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
AUX_DIR <- file.path(OUT_DIR, "analysis_snapshot")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUX_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------
required_pkgs <- c(
  "DBI",
  "duckdb",
  "data.table",
  "survival",
  "splines",
  "ggplot2",
  "scales",
  "flextable",
  "officer"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  install.packages(missing_pkgs)
}

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(data.table)
  library(survival)
  library(splines)
  library(ggplot2)
  library(scales)
  library(flextable)
  library(officer)
})

HELPER_FILE <- file.path(R_DIR, "final_analysis_helpers.R")

if (!file.exists(HELPER_FILE)) {
  stop(
    "Missing final_analysis_helpers.R at:\n",
    HELPER_FILE
  )
}

source(HELPER_FILE)

# ------------------------------------------------------------
# 2. FILE HELPERS
# ------------------------------------------------------------
first_existing <- function(paths, required = TRUE) {
  hit <- paths[file.exists(paths)]

  if (length(hit) > 0L) {
    return(hit[1])
  }

  if (required) {
    stop(
      "Could not find any of:\n",
      paste0("  - ", paths, collapse = "\n")
    )
  }

  NA_character_
}

read_required <- function(paths) {
  fread(first_existing(paths))
}

m06 <- read_required(c(
  file.path(ANALYSIS_DIR, "06_model_results.csv"),
  file.path(WORK_DIR, "06_model_results.csv")
))

m10 <- read_required(c(
  file.path(FINAL_DIR, "10_final_model_results.csv")
))

inc10 <- read_required(c(
  file.path(FINAL_DIR, "10_final_incidence_rates.csv")
))

latency <- read_required(c(
  file.path(FINAL_DIR, "10_latency_interval_results.csv")
))

competing <- read_required(c(
  file.path(FINAL_DIR, "10_competing_risk_result.csv")
))

m12 <- read_required(c(
  file.path(SIGNAL_DIR, "12_signal_model_results.csv")
))

inc12 <- read_required(c(
  file.path(SIGNAL_DIR, "12_signal_incidence_rates.csv")
))

m15 <- read_required(c(
  file.path(METHOD_DIR, "15_methodologic_model_results.csv")
))

inc15 <- read_required(c(
  file.path(METHOD_DIR, "15_methodologic_incidence_rates.csv")
))

ipcw15 <- read_required(c(
  file.path(METHOD_DIR, "15_ipcw_results.csv")
))

m16 <- read_required(c(
  file.path(STAGE16_DIR, "16_model_results.csv")
))

inc16 <- read_required(c(
  file.path(STAGE16_DIR, "16_incidence_rates.csv")
))

fdr16 <- read_required(c(
  file.path(STAGE16_DIR, "16_expanded_exploratory_FDR.csv")
))

boundary16 <- read_required(c(
  file.path(
    STAGE16_QC_DIR,
    "16_landmark_boundary_reconciliation.csv"
  )
))

boundary_audit16 <- read_required(c(
  file.path(
    STAGE16_QC_DIR,
    "16_exact_landmark_tumor_audit.csv"
  )
))

target_trial14 <- read_required(c(
  file.path(
    METHOD_QC_DIR,
    "14_target_trial_specification.csv"
  )
))

phenotype_patient14 <- read_required(c(
  file.path(
    METHOD_QC_DIR,
    "14_phenotype_corroboration_patient_level.csv"
  )
))

flow <- read_required(c(
  file.path(QC_DIR, "02_cohort_flow.csv"),
  file.path(WORK_DIR, "02_cohort_flow.csv")
))

revbase <- read_required(c(
  file.path(QC_DIR, "05_revised_base_counts.csv"),
  file.path(WORK_DIR, "05_revised_base_counts.csv")
))

# ------------------------------------------------------------
# 3. FORMAT HELPERS
# ------------------------------------------------------------
fmt_n <- function(x) {
  ifelse(
    is.na(x),
    "\u2014",
    format(
      round(x),
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )
  )
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "\u2014",
    ifelse(
      p < 0.001,
      "<0.001",
      ifelse(
        p > 0.9,
        ">0.9",
        sprintf("%.3f", p)
      )
    )
  )
}

# q values close to the 0.05 boundary need extra displayed precision:
# a value such as 0.0498 would otherwise print 0.050 yet be flagged Yes.
# Statistical decisions always use the unrounded numeric value.
fmt_q <- function(q) {
  vapply(q, function(x) {
    if (is.na(x)) return("\u2014")
    if (x < 0.001) return("<0.001")
    if (x > 0.9) return(">0.9")
    if (abs(x - 0.05) <= 0.002) {
      # Print substantial precision for borderline results.
      return(format(x, digits = 12, trim = TRUE, scientific = FALSE))
    }
    sprintf("%.3f", x)
  }, character(1))
}

fmt_hr_ci <- function(hr, lo, hi) {
  ifelse(
    is.na(hr) | is.na(lo) | is.na(hi),
    "\u2014",
    sprintf(
      "%.2f (%.2f\u2013%.2f)",
      hr,
      lo,
      hi
    )
  )
}

fmt_rate <- function(x) {
  ifelse(
    is.na(x),
    "\u2014",
    sprintf("%.1f", x)
  )
}

fmt_pct <- function(x) {
  ifelse(
    is.na(x),
    "\u2014",
    sprintf("%.1f", x)
  )
}

fmt_mean_sd <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return("\u2014")
  }

  sprintf(
    "%.1f (%.1f)",
    mean(x),
    sd(x)
  )
}

fmt_median_iqr <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return("\u2014")
  }

  q <- quantile(
    x,
    c(0.25, 0.50, 0.75),
    na.rm = TRUE
  )

  sprintf(
    "%.1f [%.1f\u2013%.1f]",
    q[2],
    q[1],
    q[3]
  )
}

fmt_n_pct <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0L) {
    return("\u2014")
  }

  n <- sum(x == 1, na.rm = TRUE)
  p <- 100 * mean(x == 1, na.rm = TRUE)

  sprintf(
    "%s (%.1f)",
    format(
      n,
      big.mark = ",",
      trim = TRUE
    ),
    p
  )
}

safe_smd <- function(x, z, w) {
  out <- try(
    smd_numeric(
      as.numeric(x),
      z = z,
      w = w
    ),
    silent = TRUE
  )

  if (
    inherits(out, "try-error") ||
    !is.finite(out)
  ) {
    return(NA_real_)
  }

  abs(as.numeric(out))
}

pretty_covariate <- function(x) {
  x <- gsub("^race_f", "Race: ", x)
  x <- gsub("^ethnicity_f", "Ethnicity: ", x)
  x <- gsub("^marital_f", "Marital status: ", x)
  x <- gsub("^region_f", "Region: ", x)
  x <- gsub("^source_f", "Source: ", x)

  lookup <- c(
    age_at_index = "Age",
    index_year = "Index year",
    hypertension = "Hypertension",
    diabetes = "Diabetes",
    obesity = "Obesity",
    ckd = "Chronic kidney disease",
    smoking = "Smoking",
    cad = "Coronary artery disease",
    heart_failure = "Heart failure",
    cerebrovascular_disease = "Cerebrovascular disease",
    hyperlipidemia = "Hyperlipidemia",
    osteoporosis = "Osteoporosis",
    alopecia = "Alopecia",
    urinary_retention = "Urinary retention",
    uti = "Urinary tract infection",
    hematuria = "Hematuria",
    bladder_stone = "Bladder stone",
    log_enc_365 = "Encounters in prior year",
    log_enc_1095 = "Encounters in prior 3 years",
    log_amb_365 = "Ambulatory encounters in prior year",
    log_ed_365 = "ED encounters in prior year",
    log_ip_365 = "Inpatient encounters in prior year"
  )

  idx <- match(x, names(lookup))
  x[!is.na(idx)] <- unname(
    lookup[idx[!is.na(idx)]]
  )

  x
}

# ------------------------------------------------------------
# 4. WORD TABLE STYLE
# ------------------------------------------------------------
make_ft <- function(
    dt,
    body_size = 10,
    header_size = 10
) {

  ft <- flextable(dt)
  ft <- theme_booktabs(ft)
  ft <- font(
    ft,
    fontname = "Arial",
    part = "all"
  )
  ft <- fontsize(
    ft,
    size = body_size,
    part = "body"
  )
  ft <- fontsize(
    ft,
    size = header_size,
    part = "header"
  )
  ft <- bold(
    ft,
    bold = TRUE,
    part = "header"
  )
  ft <- align(
    ft,
    align = "left",
    part = "all"
  )
  ft <- valign(
    ft,
    valign = "center",
    part = "all"
  )
  ft <- padding(
    ft,
    padding.top = 2,
    padding.bottom = 2,
    padding.left = 3,
    padding.right = 3,
    part = "all"
  )
  ft <- set_table_properties(
    ft,
    layout = "autofit",
    width = 1
  )

  ft
}

title_fpar <- function(number, title) {
  fpar(
    ftext(
      paste0("Table ", number, ". "),
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 12,
        italic = TRUE
      )
    ),
    ftext(
      title,
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 12
      )
    )
  )
}

note_fpar <- function(note) {
  fpar(
    ftext(
      note,
      prop = fp_text(
        font.family = "Arial",
        font.size = 8
      )
    )
  )
}

add_table_block <- function(
    doc,
    number,
    title,
    ft,
    note = NULL
) {

  doc <- body_add_fpar(
    doc,
    title_fpar(number, title)
  )

  doc <- body_add_flextable(
    doc,
    ft
  )

  if (!is.null(note)) {
    doc <- body_add_fpar(
      doc,
      note_fpar(note)
    )
  }

  doc <- body_add_par(doc, "")

  doc
}

# ------------------------------------------------------------
# 5. CANONICAL BOUNDARY-FIXED PRIMARY COHORT
# ------------------------------------------------------------
con <- dbConnect(
  duckdb::duckdb(),
  dbdir = DB_PATH,
  read_only = TRUE
)

on.exit(
  try(
    dbDisconnect(
      con,
      shutdown = TRUE
    ),
    silent = TRUE
  ),
  add = TRUE
)

if (
  !"analysis_cohort_primary_boundary_fixed" %in%
  dbListTables(con)
) {
  stop(
    "Missing analysis_cohort_primary_boundary_fixed. ",
    "Run Stage 16 first."
  )
}

primary <- as.data.table(
  dbGetQuery(
    con,
    "SELECT * FROM analysis_cohort_primary_boundary_fixed;"
  )
)

primary[
  ,
  trt := as.integer(
    exposure == "5ARI"
  )
]

# Canonical overlap weights for Table 1 and Love plot.
cat("\nBuilding canonical boundary-fixed overlap weights...\n")

wt_primary <- estimate_overlap_weights(
  primary,
  label = "Final canonical primary cohort"
)

d <- wt_primary$data
primary_balance <- wt_primary$balance

z <- d$trt
wu <- rep(1, nrow(d))
ww <- d$ow

# ------------------------------------------------------------
# 6. FINAL BOUNDARY-ALIGNED SENSITIVITY REFITS
# ------------------------------------------------------------
#
# These are already-defined analyses. They are refit only because the
# one exact-landmark D32.9 patient was present in the older Stage-10/15
# denominators.
#
primary[
  ,
  index_date := as.Date(index_date)
]

final_specific <- run_ato_generic(
  copy(primary),
  "Final boundary-fixed high-specificity intracranial meningioma",
  "meningioma_event_specific_boundary",
  "followup_days_specific_boundary"
)

final_2012 <- run_ato_generic(
  copy(primary)[
    index_date >= as.Date("2012-01-01")
  ],
  "Final boundary-fixed primary cohort: index 2012 onward",
  "meningioma_event_strict_boundary",
  "followup_days_strict_boundary"
)

final_icd10 <- run_ato_generic(
  copy(primary)[
    index_date >= as.Date("2015-10-01")
  ],
  "Final boundary-fixed primary cohort: ICD-10 era index",
  "meningioma_event_strict_boundary",
  "followup_days_strict_boundary"
)

# Boundary-align the previously defined 5-y administrative horizon.
# Risk begins at the 1-y landmark; follow-up is truncated at 5 years
# after risk-set entry. No new hypothesis is introduced.
primary_5y <- copy(primary)
five_year_days <- 5 * 365.25

primary_5y[
  ,
  event_5y_boundary := as.integer(
    meningioma_event_strict_boundary == 1L &
    followup_days_strict_boundary <= five_year_days
  )
]

primary_5y[
  ,
  followup_days_5y_boundary := pmin(
    as.numeric(followup_days_strict_boundary),
    five_year_days
  )
]

final_5y <- run_ato_generic(
  primary_5y,
  "Final boundary-fixed 5-y administrative horizon",
  "event_5y_boundary",
  "followup_days_5y_boundary"
)

final_refits <- rbindlist(
  list(
    final_specific$result,
    final_2012$result,
    final_icd10$result,
    final_5y$result
  ),
  fill = TRUE
)

final_refit_incidence <- rbindlist(
  list(
    final_specific$incidence,
    final_2012$incidence,
    final_icd10$incidence,
    final_5y$incidence
  ),
  fill = TRUE
)

fwrite(
  final_refits,
  file.path(
    AUX_DIR,
    "FINAL_boundary_aligned_sensitivity_models.csv"
  )
)

fwrite(
  final_refit_incidence,
  file.path(
    AUX_DIR,
    "FINAL_boundary_aligned_sensitivity_incidence.csv"
  )
)

# ------------------------------------------------------------
# 7. TABLE 1: BASELINE CHARACTERISTICS
# ------------------------------------------------------------
table1_rows <- list()

add_cont_row <- function(
    label,
    var,
    summary = c("mean_sd", "median_iqr")
) {

  summary <- match.arg(summary)
  x <- as.numeric(d[[var]])

  v1 <- if (summary == "mean_sd") {
    fmt_mean_sd(x[z == 1])
  } else {
    fmt_median_iqr(x[z == 1])
  }

  v0 <- if (summary == "mean_sd") {
    fmt_mean_sd(x[z == 0])
  } else {
    fmt_median_iqr(x[z == 0])
  }

  table1_rows[[length(table1_rows) + 1L]] <<-
    data.table(
      Characteristic = label,
      `5-ARI` = v1,
      `Uroselective alpha-blocker` = v0,
      `SMD before` = sprintf(
        "%.3f",
        safe_smd(x, z, wu)
      ),
      `SMD after` = sprintf(
        "%.3f",
        safe_smd(x, z, ww)
      )
    )
}

add_binary_row <- function(
    label,
    var
) {

  x <- as.integer(d[[var]])

  table1_rows[[length(table1_rows) + 1L]] <<-
    data.table(
      Characteristic = label,
      `5-ARI` = fmt_n_pct(
        x[z == 1]
      ),
      `Uroselective alpha-blocker` =
        fmt_n_pct(
          x[z == 0]
        ),
      `SMD before` = sprintf(
        "%.3f",
        safe_smd(x, z, wu)
      ),
      `SMD after` = sprintf(
        "%.3f",
        safe_smd(x, z, ww)
      )
    )
}

add_factor_block <- function(
    label,
    var,
    preferred_levels = NULL
) {

  x <- as.character(d[[var]])
  x[
    is.na(x) |
    trimws(x) == ""
  ] <- "Unknown"

  levs <- if (is.null(preferred_levels)) {
    sort(unique(x))
  } else {
    unique(
      c(
        preferred_levels[
          preferred_levels %in% unique(x)
        ],
        setdiff(
          sort(unique(x)),
          preferred_levels
        )
      )
    )
  }

  table1_rows[[length(table1_rows) + 1L]] <<-
    data.table(
      Characteristic = label,
      `5-ARI` = "",
      `Uroselective alpha-blocker` = "",
      `SMD before` = "",
      `SMD after` = ""
    )

  for (lv in levs) {

    b <- as.integer(x == lv)

    table1_rows[[length(table1_rows) + 1L]] <<-
      data.table(
        Characteristic =
          paste0("   ", lv),
        `5-ARI` =
          fmt_n_pct(
            b[z == 1]
          ),
        `Uroselective alpha-blocker` =
          fmt_n_pct(
            b[z == 0]
          ),
        `SMD before` =
          sprintf(
            "%.3f",
            safe_smd(b, z, wu)
          ),
        `SMD after` =
          sprintf(
            "%.3f",
            safe_smd(b, z, ww)
          )
      )
  }
}

add_cont_row(
  "Age, mean (SD), y",
  "age_at_index",
  "mean_sd"
)

add_cont_row(
  "Index year, median [IQR]",
  "index_year",
  "median_iqr"
)

add_factor_block(
  "Race, No. (%)",
  "race",
  c(
    "White",
    "Black or African American",
    "Asian",
    "American Indian or Alaska Native",
    "Native Hawaiian or Other Pacific Islander",
    "Other Race",
    "Unknown"
  )
)

add_factor_block(
  "Ethnicity, No. (%)",
  "ethnicity",
  c(
    "Not Hispanic or Latino",
    "Hispanic or Latino",
    "Unknown"
  )
)

for (x in list(
  c("Hypertension, No. (%)", "hypertension"),
  c("Diabetes, No. (%)", "diabetes"),
  c("Obesity, No. (%)", "obesity"),
  c("Chronic kidney disease, No. (%)", "ckd"),
  c("Smoking, No. (%)", "smoking"),
  c("Coronary artery disease, No. (%)", "cad"),
  c("Heart failure, No. (%)", "heart_failure"),
  c("Cerebrovascular disease, No. (%)", "cerebrovascular_disease"),
  c("Hyperlipidemia, No. (%)", "hyperlipidemia"),
  c("Osteoporosis, No. (%)", "osteoporosis"),
  c("Alopecia, No. (%)", "alopecia"),
  c("Urinary retention, No. (%)", "urinary_retention"),
  c("Urinary tract infection, No. (%)", "uti"),
  c("Hematuria, No. (%)", "hematuria"),
  c("Bladder stone, No. (%)", "bladder_stone")
)) {
  add_binary_row(x[1], x[2])
}

add_cont_row(
  "Encounters in prior year, median [IQR]",
  "n_encounters_prior_365",
  "median_iqr"
)

add_cont_row(
  "Ambulatory encounters in prior year, median [IQR]",
  "n_ambulatory_prior_365",
  "median_iqr"
)

add_cont_row(
  "ED encounters in prior year, median [IQR]",
  "n_ed_prior_365",
  "median_iqr"
)

add_cont_row(
  "Inpatient encounters in prior year, median [IQR]",
  "n_inpatient_prior_365",
  "median_iqr"
)

table1 <- rbindlist(
  table1_rows,
  fill = TRUE
)

setnames(
  table1,
  c(
    "5-ARI",
    "Uroselective alpha-blocker"
  ),
  c(
    paste0(
      "5-ARI\n(n=",
      fmt_n(sum(z == 1)),
      ")"
    ),
    paste0(
      "Uroselective alpha-blocker\n(n=",
      fmt_n(sum(z == 0)),
      ")"
    )
  )
)

fwrite(
  table1,
  file.path(
    TABLE_DIR,
    "Table_1_Baseline_characteristics.csv"
  )
)

# ------------------------------------------------------------
# 8. NORMALIZE/EXTRACT MODEL RESULTS
# ------------------------------------------------------------
normalize_models <- function(dt) {
  x <- copy(dt)

  if ("n_5ari" %in% names(x)) {
    setnames(
      x,
      c(
        "n_5ari",
        "n_uro_alpha",
        "events_5ari",
        "events_uro_alpha",
        "ESS_5ari",
        "ESS_uro_alpha"
      ),
      c(
        "n_treated",
        "n_control",
        "events_treated",
        "events_control",
        "ESS_treated",
        "ESS_control"
      ),
      skip_absent = TRUE
    )
  }

  x
}

nm06 <- normalize_models(m06)
nm10 <- normalize_models(m10)
nm12 <- normalize_models(m12)
nm15 <- normalize_models(m15)
nm16 <- normalize_models(m16)
nfinal <- normalize_models(final_refits)

get_rate <- function(
    incidence_dt,
    analysis_name,
    trt_value
) {

  x <- incidence_dt[
    analysis == analysis_name &
    trt == trt_value
  ]

  if (nrow(x) == 0L) {
    return(NA_real_)
  }

  x$weighted_rate_per_100k_py[1]
}

extract_standard_model <- function(
    model_dt,
    analysis_name,
    display,
    incidence_dt = NULL,
    treated_label = "5-ARI",
    control_label = "Uro-alpha"
) {

  x <- model_dt[
    analysis == analysis_name
  ]

  if (nrow(x) != 1L) {
    stop(
      "Expected one model row for: ",
      analysis_name,
      ". Found ",
      nrow(x),
      "."
    )
  }

  rate1 <- if (!is.null(incidence_dt)) {
    get_rate(
      incidence_dt,
      analysis_name,
      1
    )
  } else {
    NA_real_
  }

  rate0 <- if (!is.null(incidence_dt)) {
    get_rate(
      incidence_dt,
      analysis_name,
      0
    )
  } else {
    NA_real_
  }

  data.table(
    Analysis = display,
    N = x$n,
    Exposed = treated_label,
    Comparator = control_label,
    Events_exposed = x$events_treated,
    Events_comparator = x$events_control,
    Rate_exposed = rate1,
    Rate_comparator = rate0,
    HR = x$weighted_HR,
    CI_low = x$weighted_CI_low,
    CI_high = x$weighted_CI_high,
    p = x$weighted_p
  )
}

# ------------------------------------------------------------
# 9. TABLE 2: PRIMARY + ROBUSTNESS
# ------------------------------------------------------------
t2_raw <- rbindlist(
  list(
    extract_standard_model(
      nm16,
      "Boundary-fixed primary 1-y landmark",
      "Primary 1-y landmark",
      inc16
    ),

    extract_standard_model(
      nm15,
      "Target-trial aligned: day-0 ITT-like initiation effect, zero lag",
      "Target-trial aligned, day-0 ITT-like",
      inc15
    ),

    extract_standard_model(
      nm16,
      "Boundary-fixed target-trial day-0 assignment, 1-y landmark",
      "Target-trial aligned, 1-y landmark",
      inc16
    ),

    extract_standard_model(
      nfinal,
      "Final boundary-fixed high-specificity intracranial meningioma",
      "High-specificity D32.0 / ICD-9 225.2",
      final_refit_incidence
    ),

    extract_standard_model(
      nm06,
      "Sensitivity: 2-y lag",
      "2-y landmark"
    ),

    extract_standard_model(
      nfinal,
      "Final boundary-fixed 5-y administrative horizon",
      "5-y administrative follow-up",
      final_refit_incidence
    ),

    extract_standard_model(
      nfinal,
      "Final boundary-fixed primary cohort: index 2012 onward",
      "Index year 2012 onward",
      final_refit_incidence
    ),

    extract_standard_model(
      nfinal,
      "Final boundary-fixed primary cohort: ICD-10 era index",
      "ICD-10-era initiation",
      final_refit_incidence
    )
  ),
  fill = TRUE
)

# Fine-Gray competing-risk analysis.
if (nrow(competing) >= 1L) {

  t2_raw <- rbind(
    t2_raw,
    data.table(
      Analysis = "Fine-Gray model; death as competing event",
      N = competing$n[1],
      Exposed = "5-ARI",
      Comparator = "Uro-alpha",
      Events_exposed = NA_real_,
      Events_comparator = NA_real_,
      Rate_exposed = NA_real_,
      Rate_comparator = NA_real_,
      HR = competing$subdistribution_HR[1],
      CI_low = competing$CI_low[1],
      CI_high = competing$CI_high[1],
      p = competing$p[1]
    ),
    fill = TRUE
  )
}

# IPCW result: use the prespecified 1st/99th percentile truncation.
ipcw_row <- ipcw15[
  grepl(
    "1st/99th",
    analysis,
    fixed = TRUE
  )
]

if (nrow(ipcw_row) == 1L) {

  t2_raw <- rbind(
    t2_raw,
    data.table(
      Analysis = "EHR-continuity IPCW; 1st/99th percentile truncation",
      N = ipcw_row$n_patients,
      Exposed = "5-ARI",
      Comparator = "Uro-alpha",
      Events_exposed = NA_real_,
      Events_comparator = NA_real_,
      Rate_exposed = NA_real_,
      Rate_comparator = NA_real_,
      HR = ipcw_row$HR,
      CI_low = ipcw_row$CI_low,
      CI_high = ipcw_row$CI_high,
      p = ipcw_row$p
    ),
    fill = TRUE
  )
}

# Preferred manuscript order.
t2_order <- c(
  "Primary 1-y landmark",
  "Target-trial aligned, day-0 ITT-like",
  "Target-trial aligned, 1-y landmark",
  "High-specificity D32.0 / ICD-9 225.2",
  "2-y landmark",
  "5-y administrative follow-up",
  "Fine-Gray model; death as competing event",
  "EHR-continuity IPCW; 1st/99th percentile truncation",
  "Index year 2012 onward",
  "ICD-10-era initiation"
)

t2_raw[
  ,
  order_id := match(
    Analysis,
    t2_order
  )
]

setorder(
  t2_raw,
  order_id
)

t2_raw[
  ,
  order_id := NULL
]

table2 <- t2_raw[
  ,
  .(
    Analysis,
    `N` = fmt_n(N),
    `Events, exposed/comparator` =
      ifelse(
        is.na(Events_exposed),
        "\u2014",
        paste0(
          fmt_n(Events_exposed),
          "/",
          fmt_n(Events_comparator)
        )
      ),
    `Adjusted HR (95% CI)` =
      fmt_hr_ci(
        HR,
        CI_low,
        CI_high
      ),
    `P value` = fmt_p(p)
  )
]

fwrite(
  table2,
  file.path(
    TABLE_DIR,
    "Table_2_Primary_robustness_results.csv"
  )
)

# ------------------------------------------------------------
# 10. TABLE 3: EXPLORATORY PERSISTENCE
# ------------------------------------------------------------
#
# Main table shows the exposure-duration analyses most relevant to the
# biological hypothesis. The COMPLETE 16-analysis multiplicity family is
# reported in Table S6, including dutasteride and head-to-head analyses.
#
selected_stage12 <- c(
  "Sustained 1-y + landmark utilization: all 5ARI vs uro-alpha",
  "Sustained 1-y + landmark utilization: finasteride vs uro-alpha",
  "Strict sustained 1-y + landmark utilization: all 5ARI vs uro-alpha",
  "Strict sustained 1-y + landmark utilization: finasteride vs uro-alpha",
  "Persistent 3-y + landmark utilization: all 5ARI vs uro-alpha",
  "Persistent 3-y + landmark utilization: finasteride vs uro-alpha"
)

display_stage12 <- c(
  "Sustained 1-y: all 5-ARI vs uro-alpha",
  "Sustained 1-y: finasteride vs uro-alpha",
  "Strict sustained 1-y: all 5-ARI vs uro-alpha",
  "Strict sustained 1-y: finasteride vs uro-alpha",
  "Persistent 3-y: all 5-ARI vs uro-alpha",
  "Persistent 3-y: finasteride vs uro-alpha"
)

t3_broad <- rbindlist(
  Map(
    function(a, disp) {

      treated_label <- if (
        grepl(
          "finasteride",
          a,
          ignore.case = TRUE
        )
      ) {
        "Finasteride"
      } else {
        "5-ARI"
      }

      extract_standard_model(
        nm12,
        a,
        disp,
        inc12,
        treated_label = treated_label,
        control_label = "Uro-alpha"
      )
    },
    selected_stage12,
    display_stage12
  ),
  fill = TRUE
)

selected_stage16 <- c(
  "High-specificity: strict sustained 1-y, all 5ARI vs uro-alpha",
  "High-specificity: strict sustained 1-y, finasteride vs uro-alpha",
  "High-specificity: persistent 3-y, all 5ARI vs uro-alpha",
  "High-specificity: persistent 3-y, finasteride vs uro-alpha"
)

display_stage16 <- c(
  "Strict sustained 1-y: all 5-ARI; high-specificity outcome",
  "Strict sustained 1-y: finasteride; high-specificity outcome",
  "Persistent 3-y: all 5-ARI; high-specificity outcome",
  "Persistent 3-y: finasteride; high-specificity outcome"
)

t3_specific <- rbindlist(
  Map(
    function(a, disp) {

      treated_label <- if (
        grepl(
          "finasteride",
          a,
          ignore.case = TRUE
        )
      ) {
        "Finasteride"
      } else {
        "5-ARI"
      }

      extract_standard_model(
        nm16,
        a,
        disp,
        inc16,
        treated_label = treated_label,
        control_label = "Uro-alpha"
      )
    },
    selected_stage16,
    display_stage16
  ),
  fill = TRUE
)

t3_raw <- rbindlist(
  list(
    t3_broad,
    t3_specific
  ),
  fill = TRUE
)

# Map manuscript display rows back to raw analysis names for q-values.
q_map <- rbindlist(
  list(
    data.table(
      Analysis = display_stage12,
      analysis_raw = selected_stage12
    ),
    data.table(
      Analysis = display_stage16,
      analysis_raw = selected_stage16
    )
  )
)

t3_raw <- merge(
  t3_raw,
  q_map,
  by = "Analysis",
  all.x = TRUE,
  sort = FALSE
)

t3_raw <- merge(
  t3_raw,
  fdr16[
    ,
    .(
      analysis_raw = analysis,
      q_BH = q_BH_expanded
    )
  ],
  by = "analysis_raw",
  all.x = TRUE,
  sort = FALSE
)

t3_order <- c(
  display_stage12[1:2],
  display_stage12[3:4],
  display_stage16[1:2],
  display_stage12[5:6],
  display_stage16[3:4]
)

t3_raw[
  ,
  order_id := match(
    Analysis,
    t3_order
  )
]

setorder(
  t3_raw,
  order_id
)

t3_raw[
  ,
  order_id := NULL
]

table3 <- t3_raw[
  ,
  .(
    Analysis,
    `N` = fmt_n(N),
    `Events, exposed/comparator` =
      paste0(
        fmt_n(Events_exposed),
        "/",
        fmt_n(Events_comparator)
      ),
    `Adjusted HR (95% CI)` =
      fmt_hr_ci(
        HR,
        CI_low,
        CI_high
      ),
    `P value` = fmt_p(p),
    `BH q value` = fmt_q(q_BH)
  )
]

fwrite(
  table3,
  file.path(
    TABLE_DIR,
    "Table_3_Exploratory_persistence_results.csv"
  )
)

# ------------------------------------------------------------
# 11. TABLE S1: TARGET-TRIAL SPECIFICATION
# ------------------------------------------------------------
tableS1 <- copy(target_trial14)

setnames(
  tableS1,
  c(
    "component",
    "target_trial",
    "observational_emulation"
  ),
  c(
    "Component",
    "Target trial",
    "Observational emulation"
  )
)

fwrite(
  tableS1,
  file.path(
    TABLE_DIR,
    "Table_S1_Target_trial_specification.csv"
  )
)

# ------------------------------------------------------------
# 12. TABLE S2: COVARIATE BALANCE
# ------------------------------------------------------------
tableS2 <- primary_balance[
  ,
  .(
    Covariate =
      pretty_covariate(covariate),
    `SMD before weighting` =
      sprintf(
        "%.3f",
        abs_smd_unweighted
      ),
    `SMD after overlap weighting` =
      sprintf(
        "%.3f",
        abs_smd_weighted
      )
  )
][
  order(
    -as.numeric(
      `SMD before weighting`
    )
  )
]

fwrite(
  tableS2,
  file.path(
    TABLE_DIR,
    "Table_S2_Covariate_balance.csv"
  )
)

# ------------------------------------------------------------
# 13. TABLE S3: LATENCY
# ------------------------------------------------------------
latency[
  ,
  interval_display := fcase(
    interval == "0-<1 y",
    "0\u2013<1 y",
    interval == "1-<2 y",
    "1\u2013<2 y",
    interval == "2-<5 y",
    "2\u2013<5 y",
    interval == ">=5 y",
    "\u22655 y",
    default = interval
  )
]

tableS3 <- latency[
  ,
  .(
    `Time since initiation` =
      interval_display,
    `N at interval entry` =
      fmt_n(n),
    `Events` =
      fmt_n(events),
    `Events, 5-ARI/uro-alpha` =
      paste0(
        fmt_n(events_5ari),
        "/",
        fmt_n(events_uro_alpha)
      ),
    `Adjusted HR (95% CI)` =
      fmt_hr_ci(
        HR,
        CI_low,
        CI_high
      ),
    `P value` =
      fmt_p(p)
  )
]

fwrite(
  tableS3,
  file.path(
    TABLE_DIR,
    "Table_S3_Latency_results.csv"
  )
)

# ------------------------------------------------------------
# 14. TABLE S4: INTERNAL PHENOTYPE CORROBORATION
# ------------------------------------------------------------
#
# Stage-14 patient-level corroboration used the pre-boundary primary
# cohort. Remove the exact landmark patient identified in Stage 16 so
# this table reconciles to 1,456 canonical primary events.
#
boundary_primary_ids <- unique(
  boundary_audit16[
    cohort == "Primary 1-y landmark",
    patient_id
  ]
)

phen <- copy(
  phenotype_patient14[
    !patient_id %in%
      boundary_primary_ids
  ]
)

criteria <- c(
  high_specificity =
    "First incident diagnosis included D32.0 / ICD-9 225.2",
  repeat_code_ge30d =
    "Second meningioma diagnosis \u226530 d later",
  tumor_table =
    "ICD-O 953x tumor-table corroboration within \u00b1365 d",
  cancer_registry =
    "Cancer Registry corroboration within \u00b1365 d",
  resection_or_srs =
    "Meningioma resection / cranial SRS within 180 d",
  oncology_treatment_table =
    "Oncology-treatment-table corroboration within 365 d",
  independent_corroboration =
    "At least 1 independent tumor/treatment source",
  any_corroboration =
    "At least 1 corroborating criterion",
  isolated_unspecified_single_code =
    "Isolated D32.9 single-code event without corroboration"
)

tableS4_raw <- rbindlist(
  lapply(
    names(criteria),
    function(v) {

      e5 <- phen[
        exposure == "5ARI"
      ]

      eu <- phen[
        exposure == "URO_ALPHA"
      ]

      data.table(
        Criterion = criteria[[v]],
        overall_n =
          sum(
            phen[[v]] == 1L,
            na.rm = TRUE
          ),
        overall_N =
          nrow(phen),
        overall_pct =
          100 *
          mean(
            phen[[v]] == 1L,
            na.rm = TRUE
          ),
        n_5ari =
          sum(
            e5[[v]] == 1L,
            na.rm = TRUE
          ),
        N_5ari =
          nrow(e5),
        pct_5ari =
          100 *
          mean(
            e5[[v]] == 1L,
            na.rm = TRUE
          ),
        n_uro =
          sum(
            eu[[v]] == 1L,
            na.rm = TRUE
          ),
        N_uro =
          nrow(eu),
        pct_uro =
          100 *
          mean(
            eu[[v]] == 1L,
            na.rm = TRUE
          )
      )
    }
  )
)

tableS4 <- tableS4_raw[
  ,
  .(
    Criterion,
    `Overall, No. (%)` =
      sprintf(
        "%s (%.1f)",
        fmt_n(overall_n),
        overall_pct
      ),
    `5-ARI events, No. (%)` =
      sprintf(
        "%s (%.1f)",
        fmt_n(n_5ari),
        pct_5ari
      ),
    `Uro-alpha events, No. (%)` =
      sprintf(
        "%s (%.1f)",
        fmt_n(n_uro),
        pct_uro
      )
  )
]

fwrite(
  tableS4,
  file.path(
    TABLE_DIR,
    "Table_S4_Phenotype_corroboration.csv"
  )
)

# ------------------------------------------------------------
# 15. TABLE S5: PRESPECIFIED FALSIFICATION OUTCOMES
# ------------------------------------------------------------
nco_names <- c(
  "Negative control: benign cranial nerve neoplasm (D33.3 / ICD-9 225.1)",
  "Negative control: benign spinal cord neoplasm (D33.4 / ICD-9 225.3)"
)

nco_display <- c(
  "Benign cranial nerve neoplasm (D33.3 / ICD-9 225.1)",
  "Benign spinal cord neoplasm (D33.4 / ICD-9 225.3)"
)

nco_raw <- rbindlist(
  Map(
    function(a, disp) {
      extract_standard_model(
        nm15,
        a,
        disp,
        inc15
      )
    },
    nco_names,
    nco_display
  ),
  fill = TRUE
)

tableS5 <- nco_raw[
  ,
  .(
    Outcome = Analysis,
    `N` = fmt_n(N),
    `Events, exposed/comparator` =
      paste0(
        fmt_n(Events_exposed),
        "/",
        fmt_n(Events_comparator)
      ),
    `Adjusted HR (95% CI)` =
      fmt_hr_ci(
        HR,
        CI_low,
        CI_high
      ),
    `P value` =
      fmt_p(p)
  )
]

fwrite(
  tableS5,
  file.path(
    TABLE_DIR,
    "Table_S5_Falsification_outcomes.csv"
  )
)

# ------------------------------------------------------------
# 16. TABLE S6: COMPLETE EXPLORATORY FDR FAMILY
# ------------------------------------------------------------
pretty_exploratory_name <- function(x) {
  y <- x

  y <- gsub(
    " \\+ landmark utilization",
    "",
    y
  )

  y <- gsub(
    "all 5ARI",
    "all 5-ARI",
    y
  )

  y <- gsub(
    "uro-alpha",
    "uro-alpha",
    y
  )

  y
}

# Sanity check: complete 16-model family must have a unique row per model.
if (nrow(fdr16) != 16L || anyDuplicated(fdr16$analysis) ||
    anyNA(fdr16$q_BH_expanded)) {
  stop("Exploratory FDR family invalid: expected 16 distinct nonmissing q values.")
}

tableS6 <- fdr16[
  ,
  .(
    Analysis =
      pretty_exploratory_name(
        analysis
      ),
    `Adjusted HR (95% CI)` =
      fmt_hr_ci(
        weighted_HR,
        weighted_CI_low,
        weighted_CI_high
      ),
    `Raw P value` =
      fmt_p(weighted_p),
    `BH q value` =
      fmt_q(q_BH_expanded),
    `FDR q<0.05` =
      ifelse(
        q_BH_expanded < 0.05,
        "Yes",
        "No"
      )
  )
]

fwrite(
  tableS6,
  file.path(
    TABLE_DIR,
    "Table_S6_Complete_exploratory_FDR_family.csv"
  )
)

# ------------------------------------------------------------
# 17. TABLE S7: RARE / CLINICALLY CONSEQUENTIAL ENDPOINTS
# ------------------------------------------------------------
rare_labels <- c(
  "ICD-O tumor-table meningioma",
  "Cancer-registry-source ICD-O meningioma",
  "Clinically treated meningioma: resection or cranial SRS within 180 d",
  "Oncology-treatment-table ICD-O meningioma"
)

rare_display <- c(
  "ICD-O tumor-table meningioma",
  "Cancer Registry ICD-O meningioma",
  "Meningioma followed by resection/SRS within 180 d",
  "Oncology-treatment-table ICD-O meningioma"
)

rare_raw <- rbindlist(
  Map(
    function(a, disp) {
      extract_standard_model(
        nm10,
        a,
        disp,
        inc10
      )
    },
    rare_labels,
    rare_display
  ),
  fill = TRUE
)

tableS7 <- rare_raw[
  ,
  .(
    Analysis,
    `N` =
      fmt_n(N),
    `Events` =
      fmt_n(
        Events_exposed +
        Events_comparator
      ),
    `Adjusted HR (95% CI)` =
      fmt_hr_ci(
        HR,
        CI_low,
        CI_high
      ),
    `P value` =
      fmt_p(p)
  )
]

fwrite(
  tableS7,
  file.path(
    TABLE_DIR,
    "Table_S7_Rare_hard_endpoints.csv"
  )
)

# ------------------------------------------------------------
# 18. COMBINED WORD DOCUMENT
# ------------------------------------------------------------
ft1 <- make_ft(table1)
ft2 <- make_ft(table2)
ft3 <- make_ft(table3)

ftS1 <- make_ft(
  tableS1,
  body_size = 9,
  header_size = 9
)

ftS2 <- make_ft(tableS2)
ftS3 <- make_ft(tableS3)
ftS4 <- make_ft(tableS4)
ftS5 <- make_ft(tableS5)
ftS6 <- make_ft(
  tableS6,
  body_size = 9,
  header_size = 9
)
ftS7 <- make_ft(tableS7)

factor_headers <- which(
  table1$Characteristic %in%
    c(
      "Race, No. (%)",
      "Ethnicity, No. (%)"
    )
)

if (length(factor_headers) > 0L) {
  ft1 <- bold(
    ft1,
    i = factor_headers,
    bold = TRUE,
    part = "body"
  )
}

doc <- read_docx()

doc <- add_table_block(
  doc,
  "1",
  "Baseline Characteristics of the Primary Active-Comparator Cohort",
  ft1,
  paste0(
    "The canonical primary cohort used the boundary-fixed 1-year landmark definition, ",
    "excluding meningeal tumors diagnosed on or before the landmark. ",
    "Values are No. (%) unless otherwise indicated. ",
    "SMD indicates standardized mean difference. ",
    "Table 1 summarizes untransformed utilization counts; corresponding ",
    "propensity-score balance terms in Table S2 include log-transformed utilization."
  )
)

doc <- add_table_block(
  doc,
  "2",
  "Primary and Robustness Analyses of Incident Meningioma",
  ft2,
  paste0(
    "Hazard ratios are overlap-weighted Cox proportional hazards estimates unless otherwise noted. ",
    "The Fine-Gray model treats death as a competing event. ",
    "The IPCW sensitivity combines the primary treatment overlap weights with stabilized ",
    "inverse-probability-of-censoring weights truncated at the 1st and 99th percentiles."
  )
)

doc <- add_table_block(
  doc,
  "3",
  "Exploratory Sustained-Exposure and Persistence Analyses",
  ft3,
  paste0(
    "These analyses were exploratory. Sustained/persistent landmark models incorporated ",
    "healthcare utilization during the exposure-assessment period into the propensity-score model. ",
    "BH q values are Benjamini-Hochberg false-discovery-rate adjusted across the complete ",
    "16-analysis exploratory family shown in Supplementary Table 6. ",
    "Threshold decisions use unrounded q values."
  )
)

doc <- add_table_block(
  doc,
  "S1",
  "Target-Trial Specification and Observational Emulation",
  ftS1,
  paste0(
    "The target-trial-aligned sensitivity defines treatment strictly from information ",
    "available at day 0 and does not use treatment during days 1-30 to determine assignment."
  )
)

doc <- add_table_block(
  doc,
  "S2",
  "Covariate Balance Before and After Overlap Weighting",
  ftS2,
  paste0(
    "Absolute standardized mean differences are shown for the canonical boundary-fixed primary cohort. ",
    "Overlap weighting provides exact mean balance for included propensity-score covariates ",
    "up to numerical precision. Healthcare-utilization SMDs are based on ",
    "the log-transformed model covariates and therefore can differ from ",
    "the untransformed descriptive SMDs in Table 1."
  )
)

doc <- add_table_block(
  doc,
  "S3",
  "Interval-Specific Analyses by Time Since Treatment Initiation",
  ftS3,
  paste0(
    "HR indicates overlap-weighted hazard ratio. These interval-specific ",
    "risk sets are derived from a separate latency analysis and should ",
    "not be interpreted as sequential subsets of the primary 1-y landmark cohort."
  )
)

doc <- add_table_block(
  doc,
  "S4",
  "Internal Corroboration of Diagnosis-Defined Incident Meningioma",
  ftS4,
  paste0(
    "This analysis describes corroborating evidence available within Dataworks and should not be ",
    "interpreted as a formal positive predictive value validation study. ",
    "A diagnosis of D32.9 recorded exactly at risk-set entry was treated as ",
    "prevalent and excluded from the boundary-corrected primary cohort."
  )
)

doc <- add_table_block(
  doc,
  "S5",
  "Prespecified Falsification Outcomes",
  ftS5,
  paste0(
    "These outcomes were prespecified as bias-detection analyses. ",
    "The benign cranial nerve neoplasm analysis should be interpreted cautiously because the ",
    "assumption of no possible steroid-hormone-related causal pathway may not be secure."
  )
)

doc <- body_add_break(doc)

doc <- add_table_block(
  doc,
  "S6",
  "Complete Exploratory Analysis Family With False-Discovery-Rate Adjustment",
  ftS6,
  paste0(
    "Benjamini-Hochberg q values were calculated across 12 exploratory ",
    "persistence/drug-comparison models plus 4 high-specificity persistence models. ",
    "Boundary-only reruns and primary robustness models were not treated as ",
    "new hypotheses. Significance is determined from unrounded q values."
  )
)

doc <- add_table_block(
  doc,
  "S7",
  "Rare and Clinically Consequential Meningioma Endpoints",
  ftS7,
  paste0(
    "These sparse analyses are supportive sensitivity analyses and were not used to define ",
    "the primary outcome. They were fit before the final exact-landmark boundary correction ",
    "(N=501,446 vs N=501,445); the single removed patient had an isolated D32.9 diagnosis ",
    "on the landmark and no corroborating tumor, registry, surgery/SRS, or oncology-treatment evidence."
  )
)

print(
  doc,
  target = file.path(
    TABLE_DIR,
    "Manuscript_Tables_FINAL.docx"
  )
)

# ------------------------------------------------------------
# 19. FIGURE HELPERS
# ------------------------------------------------------------
save_figure <- function(
    p,
    stem,
    width,
    height
) {

  ggsave(
    file.path(
      FIG_DIR,
      paste0(stem, ".png")
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(
      FIG_DIR,
      paste0(stem, ".tiff")
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )

  ggsave(
    file.path(
      FIG_DIR,
      paste0(stem, ".pdf")
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "white"
  )
}

theme_manuscript <- function(
    base_size = 10
) {

  theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    theme(
      axis.text =
        element_text(
          size = base_size
        ),
      axis.title =
        element_text(
          size = base_size
        ),
      legend.title =
        element_blank(),
      legend.text =
        element_text(
          size = base_size
        ),
      strip.background =
        element_blank(),
      strip.text =
        element_text(
          face = "bold",
          size = base_size
        ),
      plot.margin =
        margin(
          8,
          18,
          8,
          8
        )
    )
}

# ------------------------------------------------------------
# 20. FIGURE 1: CANONICAL COHORT FLOW
# ------------------------------------------------------------
get_flow_n <- function(step_number) {
  flow[
    step == step_number,
    total
  ][1]
}

n_initial <- get_flow_n(1)
n_calendar <- get_flow_n(4)
n_history <- get_flow_n(5)
n_bph <- get_flow_n(6)
n_after_core <- get_flow_n(10)
n_post_landmark <- sum(
  revbase$patients
)

n_final <- nrow(primary)
n_final_5ari <- primary[
  exposure == "5ARI",
  .N
]
n_final_uro <- primary[
  exposure == "URO_ALPHA",
  .N
]
n_final_events <- sum(
  primary$meningioma_event_strict_boundary
)

flow_nodes <- data.table(
  y = 7:1,
  label = c(
    paste0(
      "Initial unambiguous 5-ARI or uroselective\n",
      "alpha-blocker initiators\nN = ",
      fmt_n(n_initial)
    ),
    paste0(
      "Age 40\u201389 y and calendar eligible\nN = ",
      fmt_n(n_calendar)
    ),
    paste0(
      "\u2265365 d observable history and encounter\n",
      "in the prior year\nN = ",
      fmt_n(n_history)
    ),
    paste0(
      "BPH diagnosis near treatment initiation\nN = ",
      fmt_n(n_bph)
    ),
    paste0(
      "After exclusions for prior prostate cancer,\n",
      "prior meningeal tumor,\n",
      "early opposite-class use, and known\n",
      "finasteride 1-mg initiation\nN = ",
      fmt_n(n_after_core)
    ),
    paste0(
      "Observed beyond the 1-y landmark\nN = ",
      fmt_n(n_post_landmark)
    ),
    paste0(
      "Final boundary-fixed primary cohort\nN = ",
      fmt_n(n_final),
      "; incident meningiomas = ",
      fmt_n(n_final_events),
      "\n5-ARI: ",
      fmt_n(n_final_5ari),
      "   |   Uro-alpha: ",
      fmt_n(n_final_uro)
    )
  )
)

flow_exclusions <- data.table(
  y = c(
    6.5,
    5.5,
    4.5,
    3.5,
    2.5,
    1.5
  ),
  label = c(
    paste0(
      "Excluded: ",
      fmt_n(
        n_initial -
        n_calendar
      )
    ),
    paste0(
      "Excluded: ",
      fmt_n(
        n_calendar -
        n_history
      )
    ),
    paste0(
      "Excluded: ",
      fmt_n(
        n_history -
        n_bph
      )
    ),
    paste0(
      "Excluded: ",
      fmt_n(
        n_bph -
        n_after_core
      )
    ),
    paste0(
      "Excluded: ",
      fmt_n(
        n_after_core -
        n_post_landmark
      )
    ),
    paste0(
      "Meningeal tumor during or on 1-y landmark: ",
      fmt_n(
        n_post_landmark -
        n_final
      )
    )
  )
)

p_flow <- ggplot() +
  geom_rect(
    data = flow_nodes,
    aes(
      xmin = 0.05,
      xmax = 0.95,
      ymin = y - 0.38,
      ymax = y + 0.38
    ),
    fill = "white",
    color = "black",
    linewidth = 0.45
  ) +
  geom_text(
    data = flow_nodes,
    aes(
      x = 0.5,
      y = y,
      label = label
    ),
    family = "sans",
    size = 2.82,
    lineheight = 0.94
  ) +
  geom_segment(
    data = data.table(
      y1 = 6.64:1.64,
      y2 = 6.36:1.36
    ),
    aes(
      x = 0.5,
      xend = 0.5,
      y = y1,
      yend = y2
    ),
    linewidth = 0.4,
    arrow = grid::arrow(
      length =
        grid::unit(
          0.08,
          "inches"
        ),
      type = "closed"
    )
  ) +
  geom_text(
    data = flow_exclusions,
    aes(
      x = 0.99,
      y = y,
      label = label
    ),
    hjust = 0,
    family = "sans",
    size = 2.55
  ) +
  coord_cartesian(
    xlim = c(0, 1.68),
    ylim = c(0.55, 7.45),
    clip = "off"
  ) +
  theme_void(
    base_family = "sans"
  ) +
  theme(
    plot.margin =
      margin(
        10,
        165,
        10,
        10
      )
  )

save_figure(
  p_flow,
  "Figure_1_Cohort_flow_FINAL",
  width = 7.4,
  height = 8.8
)

# ------------------------------------------------------------
# 21. FIGURE 2: PRIMARY / ROBUSTNESS FOREST
# ------------------------------------------------------------
f2 <- copy(t2_raw)

f2[
  ,
  label := factor(
    Analysis,
    levels = rev(
      t2_order[
        t2_order %in% Analysis
      ]
    )
  )
]

f2[
  ,
  estimate_text :=
    fmt_hr_ci(
      HR,
      CI_low,
      CI_high
    )
]

p_forest_primary <- ggplot(
  f2,
  aes(
    x = HR,
    y = label,
    xmin = CI_low,
    xmax = CI_high
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_errorbarh(
    height = 0,
    linewidth = 0.55
  ) +
  geom_point(
    shape = 18,
    size = 2.4
  ) +
  geom_text(
    aes(
      x = Inf,
      label = estimate_text
    ),
    hjust = -0.10,
    size = 2.85,
    family = "sans"
  ) +
  scale_x_log10(
    limits = c(0.65, 1.60),
    breaks = c(
      0.75,
      1,
      1.25,
      1.5
    ),
    labels = c(
      "0.75",
      "1",
      "1.25",
      "1.5"
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    x = "Hazard ratio",
    y = NULL
  ) +
  theme_manuscript(9.2) +
  theme(
    axis.text.y =
      element_text(
        hjust = 1,
        size = 8.4
      ),
    plot.margin =
      margin(
        8,
        112,
        8,
        8
      )
  )

save_figure(
  p_forest_primary,
  "Figure_2_Primary_robustness_forest_FINAL",
  width = 7.8,
  height = 5.8
)

# ------------------------------------------------------------
# 22. FIGURE 3: EXPLORATORY PERSISTENCE FOREST
# ------------------------------------------------------------
#
# Numeric estimates are intentionally drawn OUTSIDE the confidence-
# interval panel so long CIs cannot overlap the HR text.
#
f3_names <- c(
  "Sustained 1-y: all 5-ARI vs uro-alpha",
  "Sustained 1-y: finasteride vs uro-alpha",
  "Strict sustained 1-y: all 5-ARI vs uro-alpha",
  "Strict sustained 1-y: finasteride vs uro-alpha",
  "Strict sustained 1-y: finasteride; high-specificity outcome",
  "Persistent 3-y: all 5-ARI vs uro-alpha",
  "Persistent 3-y: finasteride vs uro-alpha",
  "Persistent 3-y: finasteride; high-specificity outcome"
)

f3 <- copy(
  t3_raw[
    Analysis %in% f3_names
  ]
)

f3[
  ,
  group := fcase(
    grepl(
      "^Sustained 1-y:",
      Analysis
    ),
    "1-y sustained exposure",
    grepl(
      "^Strict sustained 1-y:",
      Analysis
    ),
    "Strict 1-y sustained exposure",
    grepl(
      "^Persistent 3-y:",
      Analysis
    ),
    "3-y persistent exposure",
    default = "Other"
  )
]

f3[
  ,
  short_label := Analysis
]

f3[
  ,
  short_label :=
    gsub(
      "^Sustained 1-y: ",
      "",
      short_label
    )
]

f3[
  ,
  short_label :=
    gsub(
      "^Strict sustained 1-y: ",
      "",
      short_label
    )
]

f3[
  ,
  short_label :=
    gsub(
      "^Persistent 3-y: ",
      "",
      short_label
    )
]

f3[
  ,
  short_label :=
    gsub(
      "; high-specificity outcome",
      " [high-specificity]",
      short_label
    )
]

f3[
  ,
  estimate_text :=
    fmt_hr_ci(
      HR,
      CI_low,
      CI_high
    )
]

group_order <- c(
  "1-y sustained exposure",
  "Strict 1-y sustained exposure",
  "3-y persistent exposure"
)

f3[
  ,
  group := factor(
    group,
    levels = group_order
  )
]

f3[
  ,
  phenotype := ifelse(
    grepl(
      "high-specificity",
      Analysis
    ),
    "High-specificity",
    "Primary phenotype"
  )
]

# Explicit row order within the facets.
f3_label_order <- rev(c(
  "all 5-ARI vs uro-alpha",
  "finasteride vs uro-alpha",
  "all 5-ARI vs uro-alpha",
  "finasteride vs uro-alpha",
  "finasteride [high-specificity]",
  "all 5-ARI vs uro-alpha",
  "finasteride vs uro-alpha",
  "finasteride [high-specificity]"
))

# Make labels unique across facet data without changing displayed text.
f3[
  ,
  row_key := paste(
    group,
    short_label,
    sep = " | "
  )
]

row_levels <- rev(
  unique(
    f3[
      order(
        match(group, group_order),
        match(Analysis, f3_names)
      ),
      row_key
    ]
  )
)

f3[
  ,
  row_key := factor(
    row_key,
    levels = row_levels
  )
]

p_forest_signal <- ggplot(
  f3,
  aes(
    x = HR,
    y = row_key,
    xmin = CI_low,
    xmax = CI_high,
    shape = phenotype
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_errorbarh(
    height = 0,
    linewidth = 0.52
  ) +
  geom_point(
    size = 2.3
  ) +
  geom_text(
    data = f3,
    mapping = aes(
      x = Inf,
      y = row_key,
      label = estimate_text
    ),
    inherit.aes = FALSE,
    hjust = -0.10,
    size = 2.75,
    family = "sans"
  ) +
  facet_grid(
    group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "^.* \\| ",
        "",
        x
      )
    }
  ) +
  scale_shape_manual(
    values = c(
      "Primary phenotype" = 16,
      "High-specificity" = 18
    )
  ) +
  scale_x_log10(
    limits = c(0.60, 5.6),
    breaks = c(
      0.75,
      1,
      1.5,
      2,
      3,
      5
    ),
    labels = c(
      "0.75",
      "1",
      "1.5",
      "2",
      "3",
      "5"
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    x = "Hazard ratio",
    y = NULL
  ) +
  theme_manuscript(8.8) +
  theme(
    strip.placement = "outside",
    strip.text.y.left =
      element_text(
        angle = 0,
        hjust = 1,
        face = "bold"
      ),
    axis.text.y =
      element_text(
        size = 8.3
      ),
    legend.position = "bottom",
    panel.spacing.y =
      grid::unit(
        0.45,
        "lines"
      ),
    plot.margin =
      margin(
        8,
        125,
        8,
        8
      )
  )

save_figure(
  p_forest_signal,
  "Figure_3_Exploratory_persistence_forest_FINAL",
  width = 7.8,
  height = 7.6
)

# ------------------------------------------------------------
# 23. FIGURE S1: LOVE PLOT
# ------------------------------------------------------------
love <- primary_balance[
  order(
    -abs_smd_unweighted
  )
][
  1:min(.N, 25)
]

love[
  ,
  Covariate :=
    pretty_covariate(covariate)
]

love_long <- rbind(
  love[
    ,
    .(
      Covariate,
      SMD = abs_smd_unweighted,
      Timing = "Before weighting"
    )
  ],
  love[
    ,
    .(
      Covariate,
      SMD = abs_smd_weighted,
      Timing = "After overlap weighting"
    )
  ]
)

love_long[
  ,
  Timing := factor(
    Timing,
    levels = c(
      "Before weighting",
      "After overlap weighting"
    )
  )
]

love_order <- love[
  order(
    abs_smd_unweighted
  ),
  Covariate
]

love_long[
  ,
  Covariate := factor(
    Covariate,
    levels = love_order
  )
]

p_love <- ggplot(
  love_long,
  aes(
    x = SMD,
    y = Covariate,
    shape = Timing
  )
) +
  geom_vline(
    xintercept = 0.10,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    size = 2.3
  ) +
  scale_shape_manual(
    values = c(
      "Before weighting" = 1,
      "After overlap weighting" = 16
    )
  ) +
  scale_x_continuous(
    limits = c(
      0,
      max(
        0.45,
        max(
          love_long$SMD,
          na.rm = TRUE
        ) * 1.05
      )
    ),
    breaks = seq(
      0,
      0.5,
      0.1
    )
  ) +
  labs(
    x =
      "Absolute standardized mean difference",
    y = NULL
  ) +
  theme_manuscript(9.0) +
  theme(
    legend.position = "bottom"
  )

save_figure(
  p_love,
  "Figure_S1_Love_plot_FINAL",
  width = 7.2,
  height = 6.8
)

# ------------------------------------------------------------
# 24. FINAL RESULTS SNAPSHOT
# ------------------------------------------------------------
snapshot <- rbindlist(
  list(
    t2_raw[
      ,
      .(
        Category = "Primary / robustness",
        Analysis,
        N,
        Events =
          Events_exposed +
          Events_comparator,
        HR,
        CI_low,
        CI_high,
        p,
        q = NA_real_
      )
    ],

    t3_raw[
      ,
      .(
        Category = "Exploratory persistence",
        Analysis,
        N,
        Events =
          Events_exposed +
          Events_comparator,
        HR,
        CI_low,
        CI_high,
        p,
        q = q_BH
      )
    ]
  ),
  fill = TRUE
)

fwrite(
  snapshot,
  file.path(
    AUX_DIR,
    "FINAL_results_snapshot.csv"
  )
)

# ------------------------------------------------------------
# 25. FIGURE LEGENDS
# ------------------------------------------------------------
captions <- c(
  "FIGURE LEGENDS",
  "",
  "Figure 1. Cohort Construction for the Primary Active-Comparator Analysis.",
  paste0(
    "The canonical primary cohort comprised men with benign prostatic hyperplasia initiating ",
    "a 5alpha-reductase inhibitor (5-ARI) or uroselective alpha-blocker. ",
    "Outcome follow-up began after the 1-year landmark; meningeal tumors diagnosed on or before ",
    "the landmark were excluded from the risk set."
  ),
  "",
  "Figure 2. Primary and Robustness Analyses of Incident Meningioma.",
  paste0(
    "Points indicate adjusted hazard ratios and horizontal lines indicate 95% confidence intervals. ",
    "The vertical dashed line indicates a hazard ratio of 1.0. ",
    "Most estimates are overlap-weighted Cox proportional hazards models. ",
    "The Fine-Gray sensitivity treats death as a competing event, and the EHR-continuity sensitivity ",
    "combines treatment overlap weights with inverse-probability-of-censoring weights."
  ),
  "",
  "Figure 3. Exploratory Sustained-Exposure and Persistence Analyses.",
  paste0(
    "Exploratory landmark analyses evaluated repeated or persistent medication documentation ",
    "and incorporated healthcare utilization during the exposure-assessment period into the ",
    "propensity-score model. Diamonds denote analyses using the high-specificity intracranial ",
    "meningioma phenotype (ICD-10-CM D32.0 or ICD-9-CM 225.2). ",
    "Multiplicity-adjusted q values are reported in Table 3 and Supplementary Table 6."
  ),
  "",
  "Supplementary Figure 1. Covariate Balance Before and After Overlap Weighting.",
  paste0(
    "Shown are the 25 covariate terms with the largest absolute standardized mean differences ",
    "before weighting in the canonical boundary-fixed primary cohort. ",
    "The dashed line indicates an absolute standardized mean difference of 0.10."
  )
)

writeLines(
  captions,
  file.path(
    OUT_DIR,
    "MANUSCRIPT_FIGURE_LEGENDS.txt"
  )
)

# ------------------------------------------------------------
# 26. OUTPUT INDEX
# ------------------------------------------------------------
index <- c(
  "FINAL MANUSCRIPT OUTPUTS",
  "",
  "MAIN TABLES",
  "Table_1_Baseline_characteristics.csv",
  "Table_2_Primary_robustness_results.csv",
  "Table_3_Exploratory_persistence_results.csv",
  "",
  "SUPPLEMENTARY TABLES",
  "Table_S1_Target_trial_specification.csv",
  "Table_S2_Covariate_balance.csv",
  "Table_S3_Latency_results.csv",
  "Table_S4_Phenotype_corroboration.csv",
  "Table_S5_Falsification_outcomes.csv",
  "Table_S6_Complete_exploratory_FDR_family.csv",
  "Table_S7_Rare_hard_endpoints.csv",
  "",
  "WORD",
  "Manuscript_Tables_FINAL.docx",
  "",
  "FIGURES (PNG, TIFF, PDF)",
  "Figure_1_Cohort_flow_FINAL",
  "Figure_2_Primary_robustness_forest_FINAL",
  "Figure_3_Exploratory_persistence_forest_FINAL",
  "Figure_S1_Love_plot_FINAL",
  "",
  "ANALYSIS SNAPSHOT",
  "FINAL_results_snapshot.csv",
  "FINAL_boundary_aligned_sensitivity_models.csv",
  "FINAL_boundary_aligned_sensitivity_incidence.csv",
  "",
  "CAPTIONS",
  "MANUSCRIPT_FIGURE_LEGENDS.txt"
)

writeLines(
  index,
  file.path(
    OUT_DIR,
    "OUTPUT_INDEX.txt"
  )
)

cat("\n============================================================\n")
cat("FINAL MANUSCRIPT OUTPUTS COMPLETE\n")
cat("============================================================\n")
cat(
  "Output directory:\n",
  OUT_DIR,
  "\n\n",
  sep = ""
)
cat(
  paste(
    index,
    collapse = "\n"
  )
)
cat("\n")
